package rejfree.scalings;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;

import rejfree.StanUtils;
import rejfree.StanUtils.StanExecution;
import rejfree.models.normal.CompareStanRFOnNormalModel;
import rejfree.models.normal.NormalChain;
import rejfree.models.normal.NormalChainOptions;
import rejfree.models.normal.NormalChain.NormalChainModel;
import bayonet.coda.EffectiveSize;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class EstimateESSByDimensionality2 implements Runnable
{
  @Option
  public int minDim = 8;
  
  @Option
  public int maxDim = 512;

  @OptionSet(name = "model")
  public NormalChainOptions options = new NormalChainOptions();
  
  @Option
  public SamplingMethod method = SamplingMethod.BPS;
  
  @OptionSet(name = "stan")
  public StanUtils.StanOptions stanOptions = new StanUtils.StanOptions();
  
  public static enum SamplingMethod
  {
    STAN_OPTIMAL {
      @Override
      public SamplingMethodImplementation newInstance(
          EstimateESSByDimensionality2 instance)
      {
        return instance.new StanSampler(true);
      }
    },
    BPS {
      @Override
      public SamplingMethodImplementation newInstance(
          EstimateESSByDimensionality2 instance)
      {
        throw new RuntimeException();
      }
    };
    public abstract SamplingMethodImplementation newInstance(EstimateESSByDimensionality2 instance);
  }
  
  public class StanSampler implements SamplingMethodImplementation
  {
    final boolean useOptimalSettings;
    boolean ran = false;
    
    StanSampler(boolean useOpt) { this.useOptimalSettings = useOpt; }
    
    @Override
    public Collection<List<Double>> computeSamples()
    {
      ran = true;
      final int dim = chain.dim();
      if (useOptimalSettings)
      {
        
        double epsilon =
            Math.pow(2,   -5.0/4.0) *  // to have d=2 corresponding to epsilon = 1/2
            Math.pow(dim, -1.0/4.0); // from Radford Neal's HMC tutorial asymptotics
        int l = (int) (5.0 * 1.0 / epsilon);
        stanOptions.nStanWarmUps = 0;
        stanOptions.useNuts = false;
        stanOptions.useDiagMetric = false;
        stanOptions.stepSize = epsilon;
        stanOptions.stepSizeJitter = 0.1;
        stanOptions.intTime = epsilon * l;
      }
      
      StanExecution stanExec = new StanExecution(chain.stanModel(), stanOptions);
      stanExec.addInit(CompareStanRFOnNormalModel.VAR_NAME, exactSample);
      stanExec.run();
      Map<String, List<Double>> stanOutput = stanExec.parsedStanOutput();
      
      Collection<List<Double>> result = new ArrayList<List<Double>>();
      for (int d = 0; d < dim; d++)
        result.add(stanOutput.get(CompareStanRFOnNormalModel.stanVarName(d)));
      return result;
    }

    @Override
    public double computeCost_nLocalGradientEvals()
    {
      if (!ran)
        throw new RuntimeException();
      
      if (useOptimalSettings)
      {
        int dim = chain.dim();
        return stanOptions.nStanIters * (stanOptions.intTime / stanOptions.stepSize) * dim;
      }
      else
      {
        throw new RuntimeException(); // see page 12 of the CMD_STAN manual
      }
    }
    
  }
  
  public static interface SamplingMethodImplementation 
  {
    public Collection<List<Double>> computeSamples();
    
    public double computeCost_nLocalGradientEvals();
  }
  
  @SuppressWarnings("unused") // remove once BPS is implemented
  private NormalChainModel modelSpec = null;
  private NormalChain chain = null;
  private DoubleMatrix exactSample = null;

  @Override
  public void run()
  {
    OutputManager out = Results.getGlobalOutputManager();
    for (int dim = minDim; dim <= maxDim; dim *= 2)
    {
      options.nPairs = dim - 1;
      chain = new NormalChain(options);
      exactSample = chain.exactSample();
      modelSpec = chain.new NormalChainModel(exactSample.data);
      SamplingMethodImplementation sampler = method.newInstance(this);
      Collection<List<Double>> samples = sampler.computeSamples();
      SummaryStatistics essStats = new SummaryStatistics();
      SummaryStatistics relativeMeanEstimateErrorStats = new SummaryStatistics();
      SummaryStatistics relativeSDEstimateErrorStats = new SummaryStatistics();
      int dimIndex = 0;
      for (List<Double> univariateSamples : samples)
      {
        final double ess = EffectiveSize.effectiveSize(univariateSamples);
        essStats.addValue( ess );
        
        final double trueMean = 0.0;
        final double trueSD = Math.sqrt(chain.covarMatrix.get(dimIndex, dimIndex));
        SummaryStatistics univariateSummary = new SummaryStatistics();
        for (double univariateSample : univariateSamples)
          univariateSummary.addValue(univariateSample);
        final double estimatedMean = univariateSummary.getMean();
        final double estimatedSD = univariateSummary.getStandardDeviation();
        final double relativeMeanEstimateError = Math.abs(trueMean - estimatedMean) / trueSD;
        final double relativeSDEstimateError = Math.abs(trueSD - estimatedSD) / trueSD;
        relativeMeanEstimateErrorStats.addValue( relativeMeanEstimateError );
        relativeSDEstimateErrorStats.addValue( relativeSDEstimateError );

        dimIndex++;
      }
      final double nLocalGradEvals = sampler.computeCost_nLocalGradientEvals();
      out.printWrite("statisticsByDim", 
          "dim", dim, 
          "avgESS", essStats.getMean(), 
          "minESS", essStats.getMin(), 
          "maxESS", essStats.getMax(), 
          "nLocalGradEvals", nLocalGradEvals,
          "avgRelativeMeanEstimateError", relativeMeanEstimateErrorStats.getMean(),
          "avgRelativeSDEstimateErrorStats", relativeSDEstimateErrorStats.getMean());
    }
    out.close();
  }
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new EstimateESSByDimensionality2());
  }
}
