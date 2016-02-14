package rejfree.scalings;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.scalings.EstimateESSByDimensionality.ModelSpec;
import bayonet.coda.EffectiveSize;
import bayonet.rplot.PlotLine;
import blang.variables.RealVariable;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class OptimizeRefreshRate implements Runnable
{
  @Option
  public double samplingInterval = 1;
  
  @Option
  public int dim = 100;

  @Option
  public Random initRandom = new Random(1);

//  @OptionSet(name = "rf")
  public LocalRFRunnerOptions options = new LocalRFRunnerOptions();

  private ModelSpec spec;
  
  @Override
  public void run()
  {
    options.maxRunningTimeMilli = Long.MAX_VALUE;
    options.maxTrajectoryLength = 100000;
    options.maxSteps = Integer.MAX_VALUE;
    options.rfOptions.refreshmentMethod = RefreshmentMethod.RESTRICTED;
    
//    OutputManager out = Results.getGlobalOutputManager();
    spec = new EstimateESSByDimensionality.ModelSpec(dim, false);
        
    
    List<Double> ys = new ArrayList<>();
    List<Double> xs = new ArrayList<>();
    for (double refreshRate = 1.000; refreshRate < 1000.002; refreshRate *= 2)
    {
      xs.add(refreshRate);
      double essPerCompute = essPerComputeTime(refreshRate, 2).getMean();
      ys.add(essPerCompute);
      System.out.println(refreshRate + "," + essPerCompute);
//      System.out.println("repetition=" + essPerComputeTime(refreshRate, 20).getMean());
    }
    PlotLine.from(xs, ys).toPDF(Results.getFileInResultFolder("plot-" + dim + "d.pdf"));
  }
  
  public SummaryStatistics essPerComputeTime(double refreshRate, int nRepeats)
  {
    SummaryStatistics results = new SummaryStatistics();
    for (int i = 0; i < nRepeats; i++)
    {
      spec.initFromStatio(initRandom);
      options.rfOptions.refreshRate = refreshRate;
      LocalRFRunner rf = new LocalRFRunner(options);
      RealVariable monitored = spec.variables.get(0);
      rf.init(spec);
      rf.addSaveRaysProcessor(Collections.singleton(monitored));
      rf.run();
      List<Double> convertToSample = rf.saveRaysProcessor.convertToSample(monitored, samplingInterval);
      double ess = EffectiveSize.effectiveSize(convertToSample);
      double ess2 = EffectiveSize2.effectiveSize(convertToSample);
      System.out.println("" + ess +  " vs " + ess2);
      double time = rf.sampler.getNCollidedVariables() + rf.sampler.getNRefreshedVariables();  // (1000.0*ess/timems);
      double essPerTime = ess/time;
      results.addValue(essPerTime);
    }
    return results;
//    out.printWrite("essPerSec", "essPerTime", essPerTime, "ess", ess, "time", time);
  }
  

    
  
//  public LocalRFRunner optimizeRefreshRate(double initValue, LocalRFRunnerOptions otherOptions, Object model)
//  {
//    UnivariateFunction f = new UnivariateFunction() {
//
//      @Override
//      public double value(double x)
//      {
//        // TODO Auto-generated method stub
//        return 0;
//      }
//      
//    };
//    UnivariateOptimizer maximizer = new BrentOptimizer(0, 100);
//    maximizer.optimize(new UnivariateObjectiveFunction(f), GoalType.MAXIMIZE, new SearchInterval(0, 1, 100));
//    
//  }

  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new OptimizeRefreshRate());
  }
  

}
