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
  public int minModelIndex = 3;
  
  @Option
  public int maxModelIndex = 9;

  @Option
  public ModelType modelType = ModelType.ISOTROPIC_NORMAL;
  
  @Option
  public SamplingMethod method = SamplingMethod.BPS_LOCAL;
  
  @OptionSet(name = "model")
  public NormalChainOptions options = new NormalChainOptions();
  
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
  
  @Option
  public int nRepeats = 1;
  
  @Option(gloss = "For HMC, deviation from optimal epsilon scaling to investigate sensitivity "
      + "(NB: to make it negative, use '-- -0.2' for example)")
  public double perturbation = 0.0;
  
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
    BPS_LOCAL {
      @Override
      public SamplingMethodImplementation newInstance(
          EstimateESSByDimensionality2 instance)
      {
        return instance.new BPSampler(true);
      }
    },
    BPS_GLOBAL {
      @Override
      public SamplingMethodImplementation newInstance(
          EstimateESSByDimensionality2 instance)
      {
        return instance.new BPSampler(false);
      }
    };
    public abstract SamplingMethodImplementation newInstance(EstimateESSByDimensionality2 instance);
  }
  
  public class BPSampler implements SamplingMethodImplementation
  {
    LocalRFRunner rf;
    final boolean isLocal;
    
    public BPSampler(boolean isLocal)
    {
      if (!isLocal)
        throw new RuntimeException("GLOBAL not yet implemented (see TODO below)");
      this.isLocal = isLocal;
    }

    @Override
    public void compute()
    {
      LocalRFRunnerOptions options = new LocalRFRunnerOptions();
      options.maxRunningTimeMilli = Long.MAX_VALUE;
      options.maxSteps = Integer.MAX_VALUE;
      options.maxTrajectoryLength = bpsTrajLength;
      options.rfOptions.refreshmentMethod = isLocal ? RefreshmentMethod.LOCAL : RefreshmentMethod.GLOBAL;
      options.rfOptions.refreshRate = 1.0;
      options.silent = true;
      options.samplingRandom = mainRandom;
      
      rf = new LocalRFRunner(options);
      rf.init(model.getModelSpec()); // TODO: this needs to be make into global using some kind of adapter, if !isLocal
      rf.addMomentRayProcessor();
      rf.run();
    }

    @Override
    public List<Double> estimates(int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      
      List<Double> result = new ArrayList<>();
      for (RealVariable var : rf.model.getLatentVariables(RealVariable.class))
        result.add( power == 1 ? rf.momentRayProcessor.getMeanEstimate(var) : rf.momentRayProcessor.getSquaredVariableEstimate(var) );
      return result;
    }

    @Override
    public double nLocalGradientEvals()
    {
      return rf.sampler.getNCollisions() * (isLocal ? model.meanDegree() * Math.log(model.dim()) : model.dim()) + rf.sampler.getNRefreshments() * (isLocal ? model.meanDegree() : model.dim());
    }
    
  }
  
  public class HMCSampler implements SamplingMethodImplementation
  {
    int l;
    List<SummaryStatistics> sampleStatistics = new ArrayList<>();

    @Override
    public void compute()
    {
      if (options.offDiag != 0.0 || options.diag != 1.0)
        throw new RuntimeException();
      
      if (modelType != ModelType.ISOTROPIC_NORMAL)
        throw new RuntimeException();
      
      IsotropicNormalHMCEnergy target = new IsotropicNormalHMCEnergy(); 
      double epsilon =
          Math.pow(2,   -5.0/4.0) *  // to have d=2 corresponding to epsilon = 1/2
          Math.pow(model.dim(), -1.0/4.0 + perturbation); // from Radford Neal's HMC tutorial asymptotics
      l = (int) (5.0 * 1.0 / epsilon);
      DoubleMatrix sample = model.sampleExact();
      for (int i = 0; i < model.dim(); i++)
        sampleStatistics.add(new SummaryStatistics());
      SummaryStatistics acceptRate = new SummaryStatistics();
      for (int i = 0 ; i < nHMCIters; i++) 
      {
        DataStruct result = hmc.HMC.doIter(mainRandom, l, epsilon, sample, target, target, randomizeHMCPathLength);
        acceptRate.addValue(result.accept ? 1 : 0);
        sample = result.next_q;
        for (int d = 0; d < model.dim(); d++)
          sampleStatistics.get(d).addValue(sample.get(d));
      }
      out.printWrite("hmc-accept-rate", "dim", model.dim(), "acceptRate", acceptRate.getMean());
    }

    @Override
    public List<Double> estimates(int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      
      List<Double> result = new ArrayList<>();
      for (int d = 0; d < model.dim(); d++)
      {
        SummaryStatistics currentStats = sampleStatistics.get(d);
        result.add( (power == 1 ? currentStats.getSum() : currentStats.getSumsq()) /currentStats.getN());
      }
      return result;
    }

    @Override
    public double nLocalGradientEvals()
    {
      return ((double) model.dim()) * nHMCIters * (randomizeHMCPathLength ? (1.0+l)/2.0 : l); 
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
      if (useOptimalSettings)
      {
        double epsilon =
            Math.pow(2,   -5.0/4.0) *  // to have d=2 corresponding to epsilon = 1/2
            Math.pow(model.dim(), -1.0/4.0 + perturbation); // from Radford Neal's HMC tutorial asymptotics
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
      
      stanExec = model.stanExecution(); //chain.stanExecution(stanOptions);  TODO: make more general
      stanExec.output = BriefFiles.createTempFile();
      stanExec.addInit(CompareStanRFOnNormalModel.VAR_NAME, model.sampleExact());
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
      for (int d = 0; d < model.dim(); d++)
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
      
      if (useOptimalSettings)
        return stanOptions.nStanIters * (stanOptions.intTime / stanOptions.stepSize) * model.dim();
      else
      {
        double totalNumberOfLeapFrogs = 0;
        List<Double> nIters = stanExec.parsedStanOutput().get("n_leapfrog__");
        if (nIters == null) // when it is set, it is not output
          totalNumberOfLeapFrogs = (stanOptions.intTime / stanOptions.stepSize) * (stanOptions.nStanIters + stanOptions.nStanWarmUps );
        else
        {
          for (double d : nIters)
            totalNumberOfLeapFrogs += d;
        }
        return totalNumberOfLeapFrogs * model.dim();
      }
    }
  }
  
  public static interface SamplingMethodImplementation 
  {
    public void compute();
    public List<Double> estimates(int power);
    public double nLocalGradientEvals();
  }
  

//  
//  int dim;
  OutputManager out;
  
  private ModelInterface model;
  
  static interface ModelInterface
  {
    public Object getModelSpec();
    public double trueValue(int cDim, int power);
    public double optimalEstimatorVariance(int cDim, int power);
    public DoubleMatrix sampleExact();
    public int dim();
    public int meanDegree();
    public StanExecution stanExecution();
  }
  
  static enum ModelType
  {
    ISOTROPIC_NORMAL {
      @Override
      public ModelInterface createModel(int index, EstimateESSByDimensionality2 instance)
      {
        return instance.new IsotropicNormalModel(1 << index); // 2^index
      }
    },
    NORMAL_CHAIN {
      @Override
      public ModelInterface createModel(int index,
          EstimateESSByDimensionality2 instance)
      {
        return instance.new NormalChainModelAdapter(1 << index);
      }
    };
    public abstract ModelInterface createModel(int index, EstimateESSByDimensionality2 instance);
  }
  
  class NormalChainModelAdapter implements ModelInterface
  {
    final int dim;
    private final NormalChainModel modelSpec;
    private final NormalChain chain;
    private final DoubleMatrix exactSample;
    
    private NormalChainModelAdapter(int dim)
    {
      this.dim = dim;
      options.nPairs = dim - 1;
      chain = new NormalChain(options);
      exactSample = chain.exactSample();
      modelSpec = chain.new NormalChainModel(exactSample.data);
    }

    @Override
    public Object getModelSpec()
    {
      return modelSpec;
    }

    public double trueValue(int cDim, int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      return power == 1 ? 0.0 : modelVariance(cDim);
    }

    public double optimalEstimatorVariance(int cDim, int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      
      return modelVariance(cDim) * power;
    }
    
    private double modelVariance(int cDim)
    {
      return chain.covarMatrix.get(cDim, cDim);
    }

    @Override
    public DoubleMatrix sampleExact()
    {
      return exactSample;
    }

    @Override
    public int dim() { return dim; }

    @Override
    public int meanDegree()
    {
      return 3;
    }

    @Override
    public StanExecution stanExecution()
    {
      return chain.stanExecution(stanOptions);
    }
    
  }

  class IsotropicNormalModel implements ModelInterface
  {
    final int dim;
    
    private IsotropicNormalModel(int dim)
    {
      this.dim = dim;
    }

    @Override
    public Object getModelSpec()
    {
      return new EstimateESSByDimensionality.ModelSpec(dim, true);
    }

    public double trueValue(int cDim, int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      return power == 1 ? 0.0 : modelVariance(cDim);
    }

    public double optimalEstimatorVariance(int cDim, int power)
    {
      if (power != 1 && power != 2)
        throw new RuntimeException();
      
      return modelVariance(cDim) * power;
    }
    
    private double modelVariance(int cDim)
    {
      return 1.0;
    }

    @Override
    public DoubleMatrix sampleExact()
    {
      DoubleMatrix exactSample = new DoubleMatrix(dim);
      for (int i = 0; i < dim; i++)
        exactSample.put(i, mainRandom.nextGaussian());
      return exactSample;
    }

    @Override
    public int dim() { return dim; }

    @Override
    public int meanDegree()
    {
      return 1;
    }

    @Override
    public StanExecution stanExecution()
    {
      throw new RuntimeException("Not implemented");
    }
  }

  @Override
  public void run()
  {
    options.useLocal = method == SamplingMethod.BPS_LOCAL;
      
    out = Results.getGlobalOutputManager();
    options.random = mainRandom;
    int minDim = -1;
    for (int repeat = 0; repeat < nRepeats; repeat++)
    {
      System.out.println("repeat " + repeat + "/" + nRepeats);
      for (int modelIndex = minModelIndex; modelIndex <= maxModelIndex; modelIndex++)  //dim = minDim; dim <= maxDim; dim *= 2)
      {
        model = modelType.createModel(modelIndex, this);
        if (minDim == -1)
          minDim = model.dim();
        SamplingMethodImplementation sampler = method.newInstance(this);
        sampler.compute();
        
        for (int power = 1; power <= 2; power++)
        {
          List<Double> estimates = sampler.estimates(power);
          double nLocalGradEvals = sampler.nLocalGradientEvals();
          for (int cDim = 0; cDim < minDim; cDim++)
          {
            out.write("results", 
                "repeat", repeat,
                "nDim", model.dim(),
                "power", power,
                "cDim", cDim,
                "estimate", estimates.get(cDim),
                "trueValue", model.trueValue(cDim, power),
                "nLocalGradEvals", nLocalGradEvals,
                "optimalEstimatorVariance", model.optimalEstimatorVariance(cDim, power));
          }
        }
      }
    }
    out.close();
  }
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new EstimateESSByDimensionality2());
  }
}
