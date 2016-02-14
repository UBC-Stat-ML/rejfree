package rejfree.scalings;

import java.util.List;
import java.util.Map;

import org.jblas.DoubleMatrix;

import rejfree.StanUtils;
import rejfree.StanUtils.StanExecution;
import rejfree.StanUtils.StanOptions;
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

  @Option 
  public NormalChainOptions options = new NormalChainOptions();
  
  @Option
  public SamplingMethod method = SamplingMethod.BPS;
  
  @OptionSet(name = "stan")
  public StanUtils.StanOptions stanOptions = new StanUtils.StanOptions();
  
  private static int dimToInspect = 1;
  
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
    public List<Double> computeSamples()
    {
      ran = true;
      if (useOptimalSettings)
      {
        int dim = chain.dim();
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
      Map<String, List<Double>> stanOutput = stanExec.variableSamplesFromStanOutput();
      return stanOutput.get(CompareStanRFOnNormalModel.stanVarName(dimToInspect));
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
    public List<Double> computeSamples();
    
    public double computeCost_nLocalGradientEvals();
  }
  
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
      List<Double> samples = sampler.computeSamples();
      final double ess = EffectiveSize.effectiveSize(samples);
      out.printWrite("essByDim", "dim", dim, "ess", ess);
    }
    out.close();
  }
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new EstimateESSByDimensionality2());
  }
}
