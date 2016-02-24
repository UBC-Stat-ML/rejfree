package rejfree.scalings;

import hmc.DataStruct;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;

import rejfree.StanUtils;
import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.StanUtils.StanExecution;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.models.normal.CompareStanRFOnNormalModel;
import rejfree.models.normal.NormalChain;
import rejfree.models.normal.NormalChain.NormalChainModel;
import rejfree.models.normal.NormalChainOptions;
import rejfree.scalings.EstimateESSByDimensionality.IsotropicNormalHMCEnergy;
import blang.variables.RealVariable;
import briefj.BriefFiles;
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
  
  @Option
  public int nHMCIters = 1000;
  
  @Option
  public boolean randomizeHMCPathLength = true;
  
  @Option
  public double bpsTrajLength = 1000.0;
  
  @Option
  public Random mainRandom = new Random(1);
  
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
    STAN {
      @Override
      public SamplingMethodImplementation newInstance(
          EstimateESSByDimensionality2 instance)
      {
        return instance.new StanSampler(false);
      }
    },
    HMC_OPTIMAL {
      @Override
      public SamplingMethodImplementation newInstance(
          EstimateESSByDimensionality2 instance)
      {
        return instance.new HMCSampler();
      }
    },
    BPS {
      @Override
      public SamplingMethodImplementation newInstance(
          EstimateESSByDimensionality2 instance)
      {
        return instance.new BPSampler();
      }
    };
    public abstract SamplingMethodImplementation newInstance(EstimateESSByDimensionality2 instance);
  }
  
  public class BPSampler implements SamplingMethodImplementation
  {
    LocalRFRunner rf;
    final int dim = chain.dim();

    @Override
    public void compute()
    {
      LocalRFRunnerOptions options = new LocalRFRunnerOptions();
      options.maxRunningTimeMilli = Long.MAX_VALUE;
      options.maxSteps = Integer.MAX_VALUE;
      options.maxTrajectoryLength = bpsTrajLength;
      options.rfOptions.refreshmentMethod = RefreshmentMethod.GLOBAL;
      options.rfOptions.refreshRate = 1.0;
      options.silent = true;
      options.samplingRandom = mainRandom;
      
      rf = new LocalRFRunner(options);
      rf.init(modelSpec);
      rf.addMomentRayProcessor();
      rf.run();
    }

    @Override
    public List<Double> estimates(int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      
      List<Double> result = new ArrayList<>();
      for (RealVariable var : modelSpec.variables)
        result.add( power == 1 ? rf.momentRayProcessor.getMeanEstimate(var) : rf.momentRayProcessor.getSquaredVariableEstimate(var) );
      return result;
    }

    @Override
    public double nLocalGradientEvals()
    {
      return rf.sampler.getNCollidedVariables() + rf.sampler.getNRefreshedVariables();
    }
    
  }
  
  public class HMCSampler implements SamplingMethodImplementation
  {
    int l;
    List<SummaryStatistics> sampleStatistics = new ArrayList<>();
    final int dim = chain.dim();

    @Override
    public void compute()
    {
      if (options.offDiag != 0.0)
        throw new RuntimeException();
      
      IsotropicNormalHMCEnergy target = new IsotropicNormalHMCEnergy();
      double epsilon =
          Math.pow(2,   -5.0/4.0) *  // to have d=2 corresponding to epsilon = 1/2
          Math.pow(dim, -1.0/4.0); // from Radford Neal's HMC tutorial asymptotics
      l = (int) (5.0 * 1.0 / epsilon);
      DoubleMatrix sample = new DoubleMatrix(dim);
      for (int i = 0; i < dim; i++)
      {
        sample.put(i, options.random.nextGaussian());
        sampleStatistics.add(new SummaryStatistics());
      }
      for (int i = 0 ; i < nHMCIters; i++) 
      {
        DataStruct result = hmc.HMC.doIter(mainRandom, l, epsilon, sample, target, target, randomizeHMCPathLength);
        sample = result.next_q;
        for (int d = 0; d < dim; d++)
          sampleStatistics.get(d).addValue(sample.get(d));
      }
    }

    @Override
    public List<Double> estimates(int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      
      List<Double> result = new ArrayList<>();
      for (int d = 0; d < dim; d++)
      {
        SummaryStatistics currentStats = sampleStatistics.get(d);
        result.add( (power == 1 ? currentStats.getSum() : currentStats.getSumsq()) /currentStats.getN());
      }
      return result;
    }

    @Override
    public double nLocalGradientEvals()
    {
      final int dim = chain.dim();
      return ((double) dim) * nHMCIters * (randomizeHMCPathLength ? l/2.0 : l); 
    }
    
  }
  
  public class StanSampler implements SamplingMethodImplementation
  {
    final boolean useOptimalSettings;
    boolean ran = false;
    StanExecution stanExec = null;
    Map<String, SummaryStatistics> stanOutSummary = null;
    
    StanSampler(boolean useOpt) { this.useOptimalSettings = useOpt; }
    
    @Override
    public void compute()
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
        stanOptions.intTime = epsilon * l;
      }
      else
        stanOptions.saveWarmUp = true; // needed to compute running time
      
      stanOptions.rand = mainRandom;
      
      stanExec = chain.stanExecution(stanOptions);
      stanExec.output = BriefFiles.createTempFile();
      stanExec.addInit(CompareStanRFOnNormalModel.VAR_NAME, exactSample);
      stanExec.run();
      stanOutSummary  = stanExec.stanOutputToSummaryStatistics();
    }
    
    @Override
    public List<Double> estimates(int power)
    {
      if (!ran)
        throw new RuntimeException();
      
      if (power != 1 && power != 2)
        throw new RuntimeException();
      List<Double> result = new ArrayList<>();
      for (int d = 0; d < chain.dim(); d++)
      {
        SummaryStatistics currentStats = stanOutSummary.get(CompareStanRFOnNormalModel.stanVarName(d));
        result.add( (power == 1 ? currentStats.getSum() : currentStats.getSumsq()) /currentStats.getN());
      }
      return result;
    }

    @Override
    public double nLocalGradientEvals()
    {
      if (!ran)
        throw new RuntimeException();
      
      int dim = chain.dim();
      if (useOptimalSettings)
        return stanOptions.nStanIters * (stanOptions.intTime / stanOptions.stepSize) * dim;
      else
      {
        double totalNumberOfLeapFrogs = 0;
        for (double d : stanExec.parsedStanOutput().get("n_leapfrog__"))
          totalNumberOfLeapFrogs += d;
        return totalNumberOfLeapFrogs * dim;
      }
    }
  }
  
  public static interface SamplingMethodImplementation 
  {
    public void compute();
    public List<Double> estimates(int power);
    public double nLocalGradientEvals();
  }
  
  private NormalChainModel modelSpec = null;
  private NormalChain chain = null;
  private DoubleMatrix exactSample = null;

  @Override
  public void run()
  {
    OutputManager out = Results.getGlobalOutputManager();
    options.random = mainRandom;
    for (int nDim = minDim; nDim <= maxDim; nDim *= 2)
    {
      options.nPairs = nDim - 1;
      chain = new NormalChain(options);
      exactSample = chain.exactSample();
      modelSpec = chain.new NormalChainModel(exactSample.data);
      SamplingMethodImplementation sampler = method.newInstance(this);
      sampler.compute();
      
      for (int power = 1; power <= 2; power++)
      {
        List<Double> estimates = sampler.estimates(power);
        double nLocalGradEvals = sampler.nLocalGradientEvals();
        for (int cDim = 0; cDim < nDim; cDim++)
        {
          out.printWrite("results", 
              "nDim", nDim,
              "power", power,
              "cDim", cDim,
              "estimate", estimates.get(cDim),
              "trueValue", trueValue(cDim, power),
              "nLocalGradEvals", nLocalGradEvals,
              "optimalEstimatorVariance", optimalEstimatorVariance(cDim, power));
        }
      }
    }
    out.close();
  }
  
  private double trueValue(int cDim, int power)
  {
    if (power != 1 && power != 2)
      throw new RuntimeException();
    return power == 1 ? 0.0 : chain.covarMatrix.get(cDim, cDim);
  }

  private double optimalEstimatorVariance(int cDim, int power)
  {
    if (power != 1 && power != 2)
      throw new RuntimeException();
    
    double modelVariance = chain.covarMatrix.get(cDim, cDim);
    return modelVariance * power;
  }

  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new EstimateESSByDimensionality2());
  }
}
