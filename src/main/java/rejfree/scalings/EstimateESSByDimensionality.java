package rejfree.scalings;

import hmc.DataStruct;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jblas.DoubleMatrix;

import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.local.CollisionFactor;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.models.normal.IsotropicNormal;
import utils.MultiVariateObj;
import utils.Objective;
import bayonet.coda.EffectiveSize;
import bayonet.rplot.PlotHistogram;
import bayonet.rplot.PlotLine;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class EstimateESSByDimensionality implements Runnable
{
  @Option
  public boolean useSparse = false;
  
  @Option
  public double samplingInterval = 0.1;
  
  @Option
  public int minDim = 8;
  
  @Option
  public int maxDim = 512;
  
  @Option
  public int nRepeatsForSlopeEstimation = 1;
  
  @Option
  public int nRepeatsForSlopeEstimatorVariability = 100;
  
  @Option
  public int nBurnIn = 0;
  
  @Option
  public boolean randomizeHMCPathLength = false;
  
  @Option
  public Random initRandom = new Random(1);

  @Option
  public Random samplingRandom = new Random(1);
  
  @Option
  public SamplingMethod method = SamplingMethod.GLOBAL_BPS;

  @Option
  public RefreshmentMethod refreshmentMethod = RefreshmentMethod.GLOBAL;

  @Option
  public double refreshRate = 1.0;

  public static enum SamplingMethod
  {
    GLOBAL_BPS {
      @Override
      public Pair<Double, Double> measureESSAndTime(EstimateESSByDimensionality instance, int dim)
      {
        LocalRFRunnerOptions options = new LocalRFRunnerOptions();
        options.maxRunningTimeMilli = Long.MAX_VALUE;
        options.maxTrajectoryLength = 10000;
        options.maxSteps = Integer.MAX_VALUE;
        options.rfOptions.collectRate = 0.0;
        options.rfOptions.refreshmentMethod = instance.refreshmentMethod;
        options.rfOptions.refreshRate = instance.refreshRate;
        options.samplingRandom = instance.samplingRandom;
        
        ModelSpec spec = new ModelSpec(dim, instance.useSparse);
        spec.initFromStatio(instance.initRandom);
        LocalRFRunner rf = new LocalRFRunner(options);
        RealVariable monitored = spec.variables.get(0);
        rf.init(spec);
        rf.addSaveRaysProcessor(Collections.singleton(monitored));
        rf.run();
        List<Double> convertToSample = rf.saveRaysProcessor.convertToSample(monitored, instance.samplingInterval);
        double ess = EffectiveSize.effectiveSize(convertToSample);
        double time = rf.sampler.getNCollidedVariables() + rf.sampler.getNRefreshedVariables(); 
        return Pair.of(ess, time);
      }
    }, 
    HMC {

      @Override
      public Pair<Double, Double> measureESSAndTime(
          EstimateESSByDimensionality instance, int dim)
      {
        IsotropicNormalHMCEnergy target = new IsotropicNormalHMCEnergy();
        double epsilon =
            Math.pow(2,   -5.0/4.0) *  // to have d=2 corresponding to epsilon = 1/2
            Math.pow(dim, -1.0/4.0); // from Radford Neal's HMC tutorial asymptotics
        int l = (int) (5.0 * 1.0 / epsilon);
        DoubleMatrix sample = new DoubleMatrix(dim);
        for (int i = 0; i < dim; i++)
          sample.put(i, instance.initRandom.nextGaussian());
        List<Double> samples = new ArrayList<>();
        int nIters = 10000;
        for (int i = 0 ; i < nIters; i++) 
        {
          DataStruct result = hmc.HMC.doIter(instance.samplingRandom, l, epsilon, sample, target, target, instance.randomizeHMCPathLength);
          sample = result.next_q;
          samples.add(sample.get(0));
        }
        double ess = EffectiveSize.effectiveSize(samples);
        double time = dim * nIters * (instance.randomizeHMCPathLength ? l/2.0 : l); 
        return Pair.of(ess, time);
      }
      
    };
    public abstract Pair<Double,Double> measureESSAndTime(EstimateESSByDimensionality instance, int dim);
  }
  
  public static class IsotropicNormalHMCEnergy implements MultiVariateObj, Objective
  {
  
    @Override
    public double functionValue(DoubleMatrix point)
    {
      return + 0.5 * point.dot(point); 
    }
  
    @Override
    public DoubleMatrix mFunctionValue(DoubleMatrix vec)
    {
      return vec;
    }
  }
  
  @Override
  public void run()
  {
    List<Double> estimatedSlopes = new ArrayList<>();
    OutputManager out = Results.getGlobalOutputManager();
    for (int j = 0; j < nRepeatsForSlopeEstimatorVariability; j++)
    {
      List<Double> xs = new ArrayList<>();
      List<Double> ys = new ArrayList<>();
      
      SimpleRegression reg = new SimpleRegression();
      for (int dim = minDim; dim <= maxDim; dim *= 2)
      {
        SummaryStatistics stat = new SummaryStatistics();
        
        for (int i = 0; i < nRepeatsForSlopeEstimation; i++)
        {
          Pair<Double,Double> essAndTime = method.measureESSAndTime(this, dim);
          double ess = essAndTime.getLeft();
          double time =essAndTime.getRight();
          
          double essPerTime = ess/time;
          out.printWrite("essPerSec", "dim", dim, "innerRep", i, "outerRep", j, "essPerTime", essPerTime, "ess", ess, "time", time);
          if (i >= nBurnIn)
            stat.addValue(essPerTime);
        }
        double meanEss = stat.getMean();
        reg.addData(Math.log10(dim), Math.log10(meanEss));
        xs.add(Math.log10(dim));
        ys.add(Math.log10(meanEss));
      }
      out.printWrite("fit", "outerRep", j, "slope", reg.getSlope(), "intercept", reg.getIntercept(), "rSquare", reg.getRSquare());
      PlotLine.from(xs, ys).toPDF(Results.getFileInResultFolder("essPerTimeByDim-" + j + ".pdf"));
      estimatedSlopes.add(reg.getSlope());
      if (estimatedSlopes.size() > 1)
        PlotHistogram.from(estimatedSlopes).toPDF(Results.getFileInResultFolder("estimatedSlopes-estimatorSamplingDistribution.pdf"));
    }
    out.close();
  }
  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new EstimateESSByDimensionality());
  }
  
  public static class ModelSpec
  {
    public final List<RealVariable> variables = new ArrayList<RealVariable>();
    
    @DefineFactor
    public final List<CollisionFactor> factors = new ArrayList<>();
    
    public ModelSpec(int dim, boolean sparse)
    {
      for (int i = 0; i < dim; i++)
        variables.add(new RealVariable(0.0));
      if (sparse)
      {
        for (int i = 0; i < dim; i++)
          factors.add(new IsotropicNormal(Collections.singletonList(variables.get(i))));
      }
      else
        this.factors.add(new IsotropicNormal(variables));
    }
    
    public void initFromStatio(Random rand)
    {
      for (RealVariable var :  variables)
        var.setValue(rand.nextGaussian());
    }
  }
}
