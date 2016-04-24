package rejfree.models.normal;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import com.google.common.base.Stopwatch;

import rejfree.local.LocalRFRunner;
import rejfree.models.normal.BrownianBridge.LocalFactorModelSpec;
import bayonet.distributions.Gamma;
import blang.MCMCFactory.MCMCOptions;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.factors.Factor;
import blang.mcmc.Move;
import blang.mcmc.RealVariableOverRelaxedSlice;
import blang.processing.ProcessorContext;
import blang.variables.RealVariable;
import blang.variables.RealVariableProcessor;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class BjpsPrototype implements Runnable
{
  @Option
  public double priorRate = 10.0;
  
  @Option
  public double priorShape = 20.0;
  
  @Option
  public int nVariables = 100;
  
  @Option
  public double trajectoryLength = 10;
  
  @Option
  public int maxNIters = Integer.MAX_VALUE;
  
  @Option
  public double maxWallClockSec = 60.0;
  
  @Option
  public Random mainRandom = new Random(1);
  
  @Option
  public Method method = Method.NAIVE;
  enum Method { GIBBS, BJPS, NAIVE }
  
  private class ModelSpec
  {
    private BrownianBridge brownianBridge = BrownianBridge.regularlySpaced(nVariables);
    private final ProbabilityModel graphicalModel;
    
    @DefineFactor
    public final Gamma<?> prior = Gamma.on(brownianBridge.globalVariance).withRateShape(priorRate, priorShape);
    
    @DefineFactor
    public final Factor brownianLikelihood = brownianBridge.fullFactor();
    
    private ModelSpec()
    {
      this.graphicalModel = ProbabilityModel.parse(this);
    }
    
    private Move sliceSamplerOnGlobalVariance()
    {
      return graphicalModel.instantiateOperator(model.brownianBridge.globalVariance, RealVariableOverRelaxedSlice.class);
    }
    
    private List<Move> sliceSamplerOnLatentGaussianVariables()
    {
      List<Move> result = new ArrayList<>();
      LocalFactorModelSpec localFactorModelSpec = brownianBridge.localFactorModelSpec();
      ProbabilityModel gm = ProbabilityModel.parse(localFactorModelSpec);
      for (RealVariable var : brownianBridge.variables)
        result.add(gm.instantiateOperator(var, RealVariableOverRelaxedSlice.class));
      return result;
    }
    
    private LocalRFRunner bpsOnBrownianBridge()
    {
      LocalRFRunner runner = new LocalRFRunner();
      runner.init(brownianBridge.localFactorModelSpec());
      runner.options.maxSteps = Integer.MAX_VALUE;
      runner.options.maxTrajectoryLength = trajectoryLength;
      runner.options.samplingRandom = mainRandom;
      runner.options.rfOptions.collectRate = 0.0;
      runner.options.silent = true;
      return runner;
    }
  }
  
  private OutputManager out = new OutputManager();
  
  private ModelSpec model = null;
  private void init()
  {
    model = new ModelSpec();
  }
  
  @Override
  public void run()
  {
    out.setOutputFolder(Results.getResultFolder());
    {
      init();
      Estimates estimates = null;
      
      if (method == Method.GIBBS)
        estimates = gibbs();
      else if (method == Method.BJPS)
        estimates = bjps();
      else if (method == Method.NAIVE)
        estimates = naive();
      else
        throw new RuntimeException();
      
      report("mean", 0.0, estimates.mean);
      report("variance", analyticMarginalizedVariance(), estimates.variance);
    }
    out.close();
  }
  
  private void report(String statName, double truth, double estimate)
  {
    out.printWrite("results",
      "statistic", statName,
      "estimate", estimate,
      "truth", truth,
      "error", Math.abs(estimate - truth));
  }
  
  private class Estimates
  {
    double variance = Double.NaN;
    double mean = Double.NaN;
    
    final RealVariableProcessor processor = new RealVariableProcessor("param", model.brownianBridge.globalVariance);

    void processParam(int i)
    {
      if (maxNIters < Integer.MAX_VALUE) // otherwise the output will not get printed anyways and may become very large
        processor.process(new ProcessorContext(i, null, options));
    }
  }
  
  MCMCOptions options = new MCMCOptions();
  {
    options.nMCMCSweeps = maxNIters;
  }
  
  private Estimates bjps()
  {
    Estimates result = new Estimates();
    
    double sum = 0.0;
    double sumSqrs = 0.0;
    double T = 0.0;
    
    RealVariable monitored = model.brownianBridge.variables.get(monitoredVariable());
    
    while (moreIterationsNeeded())
    {
      // sample the Gaussian field
      LocalRFRunner bps = model.bpsOnBrownianBridge();
      bps.addMomentRayProcessor();
      bps.run();
      
      // sample top-level parameter
      model.sliceSamplerOnGlobalVariance().execute(mainRandom);
      
      // collect
      sum     += bps.momentRayProcessor.sum.getCount(monitored);
      sumSqrs += bps.momentRayProcessor.sumSq.getCount(monitored);
      T       += bps.momentRayProcessor.currentTime;
      result.processParam(curIter);
    }
    
    result.mean = sum / T;
    result.variance = sumSqrs / T - (result.mean * result.mean);
    
    return result;
  }
  
  int curIter = 0;
  Stopwatch watch = null;
  private boolean moreIterationsNeeded()
  {
    if (watch == null)
      watch = Stopwatch.createStarted();
    
    curIter++;
    
    if (curIter >= maxNIters)
      return false;
    
    if (watch.elapsed(TimeUnit.MILLISECONDS) > 1000.0 * maxWallClockSec)
      return false;
    
    return true;
  }

  private Estimates naive()
  {
    Estimates result = new Estimates();
    
    SummaryStatistics monitorStatistics = new SummaryStatistics();
    
    while (moreIterationsNeeded())
    {
      // sample the Gaussian field
      for (Move move : model.sliceSamplerOnLatentGaussianVariables())
        move.execute(mainRandom);
      
      // sample the top-level parameter
      model.sliceSamplerOnGlobalVariance().execute(mainRandom);
      
      // collect
      monitorStatistics.addValue(model.brownianBridge.variables.get(monitoredVariable()).getValue());
      result.processParam(curIter);
    }
    
    result.mean = monitorStatistics.getMean();
    result.variance = monitorStatistics.getVariance();
    
    return result;
  }

  private Estimates gibbs()
  {
    Estimates result = new Estimates();
    
    SummaryStatistics monitorStatistics = new SummaryStatistics();
    
    while (moreIterationsNeeded())
    {
      // sample the Gaussian field
      model.bpsOnBrownianBridge().run();
      
      // sample the top-level parameter
      model.sliceSamplerOnGlobalVariance().execute(mainRandom);
      
      // collect
      monitorStatistics.addValue(model.brownianBridge.variables.get(monitoredVariable()).getValue());
      result.processParam(curIter);
    }
    
    result.mean = monitorStatistics.getMean();
    result.variance = monitorStatistics.getVariance();
    
    return result;
  }
  
  private int monitoredVariable()
  {
    return model.brownianBridge.variables.size() / 2;
  }
  
  private double analyticMarginalizedVariance()
  {
    final double t = model.brownianBridge.ts.get(monitoredVariable());
    return t * (1.0 - t) * priorShape / priorRate;
  }
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new BjpsPrototype());
  }
}
