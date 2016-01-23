package rejfree.scalings;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.local.CollisionFactor;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.models.normal.IsotropicNormal;
import bayonet.coda.EffectiveSize;
import bayonet.rplot.PlotHistogram;
import bayonet.rplot.PlotLine;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class CompareESSLocalGlobal implements Runnable
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
  public int nRepeatsForSlopeEstimation = 2;
  
  @Option
  public int nRepeatsForSlopeEstimatorVariability = 100;
  
  @Option
  public int nBurnIn = 1;
  
  @Option
  public Random initRandom = new Random(1);

//  @OptionSet(name = "rf")
  public LocalRFRunnerOptions options = new LocalRFRunnerOptions();
  
  @Override
  public void run()
  {
    options.maxRunningTimeMilli = Long.MAX_VALUE;
    options.maxTrajectoryLength = 10000;
    options.maxSteps = Integer.MAX_VALUE;
    options.rfOptions.refreshmentMethod = RefreshmentMethod.GLOBAL;
//    options.rfOptions.refreshRate = 0.1;
    
    List<Double> estimatedSlopes = new ArrayList<>();
    OutputManager out = Results.getGlobalOutputManager();
    for (int j = 0; j < nRepeatsForSlopeEstimatorVariability; j++)
    {
      List<Double> xs = new ArrayList<>();
      List<Double> ys = new ArrayList<>();
      
      SimpleRegression reg = new SimpleRegression();
      for (int dim = minDim; dim <= maxDim; dim *= 2)
      {
//        options.rfOptions.refreshRate /= 2.0;
        SummaryStatistics stat = new SummaryStatistics();
        ModelSpec spec = new ModelSpec(dim, useSparse);
        spec.initFromStatio(initRandom);
        for (int i = 0; i < nRepeatsForSlopeEstimation; i++)
        {
          LocalRFRunner rf = new LocalRFRunner(options);
          RealVariable monitored = spec.variables.get(0);
          rf.init(spec);
          rf.addSaveRaysProcessor(Collections.singleton(monitored));
          rf.run();
          List<Double> convertToSample = rf.saveRaysProcessor.convertToSample(monitored, samplingInterval);
          double ess = EffectiveSize.effectiveSize(convertToSample);
//          long timems = rf.watch.elapsed(TimeUnit.MILLISECONDS);
          double time = rf.sampler.getNCollidedVariables() + rf.sampler.getNRefreshedVariables();  // (1000.0*ess/timems);
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
    Mains.instrumentedRun(args, new CompareESSLocalGlobal());
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
