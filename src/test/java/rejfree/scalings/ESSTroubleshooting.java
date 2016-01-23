package rejfree.scalings;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.local.CompareESSLocalGlobal;
import rejfree.local.CompareESSLocalGlobal.ModelSpec;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import bayonet.rplot.PlotLine;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class ESSTroubleshooting implements Runnable
{

  @Option
  public int order = 2;
  
  @Option
  public int nReplicates = 100;

  @Option
  public double L = 10000;

  @Option 
  public double refreshRate = 1.0;
  
  @Option
  public Random rand = new Random(1);
  
  @Option
  public int dim = 2;

  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new ESSTroubleshooting());
  }
    
  @Override
  public void run()
  {
    List<Double> xs = new ArrayList<>(), ys = new ArrayList<>();
    SimpleRegression reg = new SimpleRegression();
    OutputManager out = Results.getGlobalOutputManager();
    
    for (int d = 2; d <= 1024; d *= 2)
    {
      dim = d;
      
      
      double x = Math.log10(d);
      SummaryStatistics estimate = indirectAutocorrelationEstimator();
      
      double time = runningTimeStats.getMean();
      double ess = estimate.getMean();
      double essStdError = estimate.getStandardDeviation() / Math.sqrt(nReplicates);
      double essPerTime = ess/time;
      out.printWrite("essPerSec", "dim", dim, "essPerTime", essPerTime, "ess", ess, "essStdError", essStdError, "time", time);
      double y = Math.log10(essPerTime);
      
      reg.addData(x, y);
      xs.add(x);
      ys.add(y);
      System.out.println("slope so far: " + reg.getSlope());
      PlotLine.from(xs, ys).toPDF(Results.getFileInResultFolder("regression.pdf"));
    }
    out.printWrite("fit", "slope", reg.getSlope(), "intercept", reg.getIntercept(), "rSquare", reg.getRSquare());
    out.close();
  }
  
  SummaryStatistics runningTimeStats = null;
  public SummaryStatistics indirectAutocorrelationEstimator()
  {
    if (order > 2 || order == 0)
      throw new RuntimeException();
    
    runningTimeStats = new SummaryStatistics();
    
    SummaryStatistics result = new SummaryStatistics();
    
    double referenceVariance = order == 1 ? 1.0 : 2.0;
    
    ModelSpec spec = new CompareESSLocalGlobal.ModelSpec(dim, false);
    LocalRFRunnerOptions options = new LocalRFRunnerOptions();
    options.samplingRandom = rand;
    options.maxRunningTimeMilli = Long.MAX_VALUE;
    options.maxSteps = Integer.MAX_VALUE;
    options.maxTrajectoryLength = L;
    options.rfOptions.refreshmentMethod = RefreshmentMethod.GLOBAL;
    options.rfOptions.refreshRate = refreshRate;
    options.silent = true;
    
    for (int rep = 0; rep < nReplicates; rep++)
    {
      spec.initFromStatio(rand);
      LocalRFRunner rf = new LocalRFRunner(options);
      RealVariable var = spec.variables.get(0);
      rf.init(spec);
      rf.addMomentRayProcessor();
      rf.run();
      runningTimeStats.addValue( rf.sampler.getNCollidedVariables() + rf.sampler.getNRefreshedVariables() );
      
      double estimate = order == 1 ? rf.momentRayProcessor.getMeanEstimate(var) : rf.momentRayProcessor.getSquaredVariableEstimate(var);
      double truth = order == 1 ? 0.0 : 1.0;
      result.addValue( L * (estimate - truth) * (estimate - truth) / referenceVariance );
    }
    
    return result;
  }
}
