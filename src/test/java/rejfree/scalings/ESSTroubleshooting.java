package rejfree.scalings;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import bayonet.coda.EffectiveSize;
import bayonet.rplot.PlotLine;
import blang.variables.RealVariable;
import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.local.CompareESSLocalGlobal;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.local.CompareESSLocalGlobal.ModelSpec;



public class ESSTroubleshooting
{

  public static void main(String [] args)
  {
    
    
//    double a = 0.5;
//    double b = 0.4;
//    double c = 0.3;
//    
//    System.out.println(  a*(a-b) + a*b*(1-a)/a + (1-a)*c*a/(1-a) + (1-a)*(1-c)  );
    
    
    
    
    
    Random rand = new Random(1);
    
    // slope for # of collisions
    
    List<Double> xs = new ArrayList<>(), ys = new ArrayList<>();
    SimpleRegression reg = new SimpleRegression();
    
    
    
    for (int d = 2; d <= 1024; d *= 2)
    {
      SummaryStatistics result = new SummaryStatistics();
      
      ModelSpec spec = new CompareESSLocalGlobal.ModelSpec(d, false);
      LocalRFRunnerOptions options = new LocalRFRunnerOptions();
      options.samplingRandom = rand;
      options.maxRunningTimeMilli = Long.MAX_VALUE;
      options.maxSteps = Integer.MAX_VALUE;
      options.maxTrajectoryLength = 100;
      options.rfOptions.refreshmentMethod = RefreshmentMethod.GLOBAL;
      options.rfOptions.refreshRate = 1.0;
      options.silent = true;
      
      int nReplicates = 10;
      for (int rep = 0; rep < nReplicates; rep++)
      {
        spec.initFromStatio(rand);
        LocalRFRunner rf = new LocalRFRunner(options);
        rf.init(spec);
        rf.run();
        result.addValue(rf.sampler.getNCollisions());
      }
      System.out.println( " -> " + result.getMean() + " +/- " + result.getStandardDeviation() / Math.sqrt(nReplicates));
      double x = Math.log10(d);
      double y = Math.log10(result.getMean());
      reg.addData(x, y);
      xs.add(x);
      ys.add(y);
    }
    System.out.println("reg slope is " + reg.getSlope() + ", r2 = " + reg.getRSquare());
    PlotLine.from(xs, ys).toPDF(new File("temp.pdf"));
    
    System.out.println(indirectAutocorrelationEstimator(2,100, 100, 1000, 1.0, rand));
    
    for (int d = 2; d < 1024; d*=2)
    {
      System.out.println("d=" + d);
      System.out.println(indirectAutocorrelationEstimator(2,1000, 1000, d, 1.0, rand));
    }
    
    
    System.out.println(indirectAutocorrelationEstimator(1,100, 1000, 2, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(1,100, 1000, 2, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(1,100, 1000, 2, 1.0, rand));
    System.out.println("--");
    System.out.println(indirectAutocorrelationEstimator(1,1000, 100, 2, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(1,10000, 100, 2, 1.0, rand));
    
    System.out.println(indirectAutocorrelationEstimator(2,100, 1000, 2, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(2,1000, 1000, 2, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(2,10000, 1000, 2, 1.0, rand));
    
    
    List<Double> indeps = new ArrayList<Double>();
    for (int i = 0; i < 100; i++)
      indeps.add(rand.nextGaussian());
    System.out.println(EffectiveSize.effectiveSize(indeps));
    System.out.println(EffectiveSize2.effectiveSize(indeps));
    
    
    System.out.println(indirectAutocorrelationEstimator(2,100, 1000, 2, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(2,100, 1000, 20, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(2,100, 1000, 200, 1.0, rand));
    
    System.out.println(indirectAutocorrelationEstimator(2,100, 1000, 100, 0.1, rand));
    System.out.println(indirectAutocorrelationEstimator(2,100, 1000, 100, 1.0, rand));
    System.out.println(indirectAutocorrelationEstimator(2,100, 1000, 100, 10, rand));
  }
  
  public static double indirectAutocorrelationEstimator(int order, int nReplicates, double L, int dim, double refreshRate, Random rand)
  {
    if (order > 2 || order == 0)
      throw new RuntimeException();
    
    SummaryStatistics result = new SummaryStatistics();
    
    ModelSpec spec = new CompareESSLocalGlobal.ModelSpec(dim, false);
    LocalRFRunnerOptions options = new LocalRFRunnerOptions();
    options.samplingRandom = rand;
    options.maxRunningTimeMilli = Long.MAX_VALUE;
    options.maxSteps = Integer.MAX_VALUE;
    options.maxTrajectoryLength = L;
    options.rfOptions.refreshmentMethod = RefreshmentMethod.GLOBAL;
    options.rfOptions.refreshRate = refreshRate;
    options.silent = true;
    
//    SummaryStatistics nColls = new SummaryStatistics();
//    SummaryStatistics ess1 = new SummaryStatistics();
//    SummaryStatistics ess2 = new SummaryStatistics();
    
    for (int rep = 0; rep < nReplicates; rep++)
    {
      spec.initFromStatio(rand);
      LocalRFRunner rf = new LocalRFRunner(options);
      RealVariable var = spec.variables.get(0);
      rf.init(spec);
      rf.addMomentRayProcessor();
//      rf.addSaveRaysProcessor(Collections.singletonList(var));
      rf.run();
//      nColls.addValue(rf.sampler.getNCollisions());
      
//      List<Double> convertToSample = rf.saveRaysProcessor.convertToSample(var, 1);
//      ess1.addValue(EffectiveSize.effectiveSize(convertToSample));
//      ess2.addValue(EffectiveSize2.effectiveSize(convertToSample));
      
      double estimate = order == 1 ? rf.momentRayProcessor.getMeanEstimate(var) : rf.momentRayProcessor.getSquaredVariableEstimate(var);
      double truth = order == 1 ? 0.0 : 1.0;
      result.addValue( (estimate - truth) * (estimate - truth) );
    }
    
//    System.out.println("nColl=" + nColls.getMean());
//    System.out.println("acf1=" + (L/ess1.getMean()));
//    System.out.println("acf1=" + (L/ess2.getMean()));
    
    double referenceVariance = order == 1 ? 1.0 : 2.0;
    
    
    System.out.println( " -> " + result.getMean() + " +/- " + result.getStandardDeviation() / Math.sqrt(nReplicates));
    
    return L * result.getMean() / referenceVariance;
  }
}
