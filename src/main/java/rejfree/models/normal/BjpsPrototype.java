package rejfree.models.normal;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import rejfree.local.LocalRFRunner;
import bayonet.distributions.Gamma;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.MCMCFactory.MCMCOptions;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.factors.Factor;
import blang.mcmc.Move;
import blang.mcmc.RealVariableOverRelaxedSlice;
import blang.mcmc.RealVariablePeskunTypeMove;
import blang.processing.Processor;
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
  public int nVariables = 1;
  
  @Option
  public double trajectoryLength = 1.0;
  
  @Option
  public int nIters = 10_000;
  
  @Option
  public Random mainRandom = new Random(1);
  
  @Option
  public Method method = Method.NAIVE_MH;
  enum Method { GIBBS, BJPS, NAIVE_MH }
  
  private ModelSpec model = null;

  private void init()
  {
    model = new ModelSpec(method == Method.BJPS);
  }
  
  private class ModelSpec
  {
    private BrownianBridge brownianBridge = BrownianBridge.regularlySpaced(nVariables);
    private final ProbabilityModel graphicalModel;
    
    @DefineFactor
    public final Gamma<?> prior = Gamma.on(brownianBridge.globalVariance).withRateShape(priorRate, priorShape);
    
    @DefineFactor
    public final Factor brownianLikelihood;
    
    @DefineFactor 
    public final List<? extends Factor> localFactors;
    
    @DefineFactor
    public final Factor priorNormalizationFactor;
    
    private ModelSpec(boolean useFullSlice)
    {
      brownianLikelihood = useFullSlice ? null : brownianBridge.fullFactor();
      localFactors = useFullSlice ? brownianBridge.localCollisionFactors() : null;
      priorNormalizationFactor = useFullSlice ? brownianBridge.normalizationFactor() : null;
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
      Estimates estimates = null;
      
      if (method == Method.GIBBS)
        estimates = gibbs();
      else if (method == Method.BJPS)
        estimates = bjps();
      else if (method == Method.NAIVE_MH)
        estimates = fullSlice();
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
  
  private Estimates fullSlice()
  {
    Estimates result = new Estimates();
    SummaryStatistics monitorStatistics = new SummaryStatistics();
    
    MCMCFactory factory = new MCMCFactory();
    
    {
      /*
       * Found bug in bayonet: the sampler RealVariableOverRelaxedSlice 
       * appears not to be pi-invariant in this model. 
       * 
       * To reproduce the problem make the first line below uncommented 
       * and the second one commented. The error does not go to zero.
       */
      
  //    factory.excludeNodeMove(RealVariablePeskunTypeMove.class); 
      factory.excludeNodeMove(RealVariableOverRelaxedSlice.class);
    }
    
    factory.excludeNodeProcessor(RealVariableProcessor.class);
    factory.mcmcOptions.nMCMCSweeps = nIters;
    factory.mcmcOptions.burnIn = 0;
    factory.mcmcOptions.thinningPeriod = 1;
    factory.mcmcOptions.random = mainRandom;
    factory.addProcessor(new Processor() {
      @Override
      public void process(ProcessorContext context)
      {
        monitorStatistics.addValue(model.brownianBridge.variables.get(monitoredVariable()).getValue());
        result.processParam(context.getMcmcIteration());
      }
    });
    
    MCMCAlgorithm algo = factory.build(model.graphicalModel);
    algo.run();
    
    result.mean = monitorStatistics.getMean();
    result.variance = monitorStatistics.getVariance();
    
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
