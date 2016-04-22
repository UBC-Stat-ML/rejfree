package rejfree.models.normal;

import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import rejfree.local.LocalRFRunner;
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
  public double priorRate = 1.0;
  
  @Option
  public double priorShape = 2.0;
  
  @Option
  public int nVariables = 1000;
  
  @Option
  public double trajectoryLength = 1.0;
  
  @Option
  public int nIters = 10000;
  
  @Option
  public Random mainRandom = new Random(1);
  
  @Option
  public Method method = Method.BJPS;
  enum Method { GIBBS, BJPS }
  
  private ModelSpec model = null;

  private void init()
  {
    model = new ModelSpec();
  }
  
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
  
  OutputManager out = new OutputManager();
  
  @Override
  public void run()
  {
    out.setOutputFolder(Results.getResultFolder());
    {
      init();
      Estimates estimates = useGibbs() ? gibbs() : bjps();
      report("mean", 0.0, estimates.mean);
      report("variance", analyticMarginalizedVariance(), estimates.variance);
    }
    out.close();
  }
  
  private boolean useGibbs()
  {
    if (method == Method.GIBBS)
      return true;
    else if (method == Method.BJPS)
      return false;
    else
      throw new RuntimeException();
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
      processor.process(new ProcessorContext(i, null, options));
    }
  }
  
  MCMCOptions options = new MCMCOptions();
  {
    options.nMCMCSweeps = nIters;
  }
  
  private Estimates bjps()
  {
    Estimates result = new Estimates();
    
    double sum = 0.0;
    double sumSqrs = 0.0;
    double T = 0.0;
    
    RealVariable monitored = model.brownianBridge.variables.get(monitoredVariable());
    
    for (int i = 0; i < nIters; i++)
    {
      model.sliceSamplerOnGlobalVariance().execute(mainRandom);
      LocalRFRunner bps = model.bpsOnBrownianBridge();
      bps.addMomentRayProcessor();
      bps.run();
      // collect
      sum     += bps.momentRayProcessor.sum.getCount(monitored);
      sumSqrs += bps.momentRayProcessor.sumSq.getCount(monitored);
      T       += bps.momentRayProcessor.currentTime;
      result.processParam(i);
    }
    
    result.mean = sum / T;
    result.variance = sumSqrs / T - (result.mean * result.mean);
    
    return result;
  }
  
  private Estimates gibbs()
  {
    Estimates result = new Estimates();
    
    SummaryStatistics monitorStatistics = new SummaryStatistics();
    
    for (int i = 0; i < nIters; i++)
    {
      model.sliceSamplerOnGlobalVariance().execute(mainRandom);
      model.bpsOnBrownianBridge().run();
      // collect
      monitorStatistics.addValue(model.brownianBridge.variables.get(monitoredVariable()).getValue());
      result.processParam(i);
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
