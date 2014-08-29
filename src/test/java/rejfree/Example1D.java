package rejfree;

import java.io.File;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import bayonet.distributions.Poisson;
import bayonet.distributions.Uniform;
import bayonet.rplot.PlotHistogram;
import briefj.collections.Counter;

import com.google.common.collect.Lists;



public class Example1D
{
  private static Counter<Integer> errors;
  private static int nIters = 100;
  
  public static void main(String [] args)
  {
    errors = new Counter<Integer>();
  }
  
  public static List<Double> method1(Random rand)
  {
    List<Double> result = Lists.newArrayList();
    double meanEstimator = 0.0;
    double meanEstimatorDenom = 0.0;
    for (int i = 0; i < nIters; i++)
    {
      double u = rand.nextDouble();
      double t = Math.sqrt(-2.0 * Math.log(u));
      meanEstimator += t * Math.abs(t) / 2.0;
      meanEstimatorDenom += t;
      int pp = Poisson.generate(rand, t);
      for (int j = 0; j < pp; j++)
        result.add(Uniform.generate(rand, -t, t));
    }
    
//    System.out.println("RB:" + meanEstimator/meanEstimatorDenom);
    
    return result;
  }
  
  public static List<Double> method2(Random rand)
  {
    List<Double> result = Lists.newArrayList();
    
    for (int i = 0; i < nIters; i++)
      result.add((rand.nextGaussian()));
    
    return result;
  }
  
//  public static double sampleIndepStdNormal(Random rand)
//  {
//    while (true)
//    {
//      double u1 = rand.nextDouble();
//      double t = Math.sqrt(-2.0 * Math.log(u1));
//      double u2 = rand.nextDouble();
//      if (u2 > Math.exp(-t))
//        return Uniform.generate(rand, -t, t);
//    }
//  }
//  
//  public static void main(String [] args)
//  {
//    Random rand = new Random(1);
//    int nPoints = 1000000;
//    
//    {
//      List<Double> points = Lists.newArrayList();
//      SummaryStatistics stat = new SummaryStatistics();
//      for (int i = 0; i < nPoints; i++)
//      {
//        double current = sampleIndepStdNormal(rand);
//        stat.addValue(current);
//        points.add(current);
//      }
//      System.out.println(stat);
//      PlotHistogram.from(points).toPDF(new File("hist.pdf"));
//    }
//    {
//      List<Double> points = Lists.newArrayList();
//      SummaryStatistics stat = new SummaryStatistics();
//      for (int i = 0; i < nPoints; i++)
//      {
//        double current = rand.nextGaussian();
//        stat.addValue(current);
//        points.add(current);
//      }
//      System.out.println(stat);
//      PlotHistogram.from(points).toPDF(new File("hist2.pdf"));
//    }
//    
//    System.out.println("---");
//    
////    Random rand = new Random(1);
//    summarize(method1(rand));
//    summarize(method2(rand));
//  }
//  
//  public static void summarize(List<Double> samples)
//  {
//    SummaryStatistics stats = new SummaryStatistics();
//    for (double d : samples)
//      stats.addValue(d);
//    System.out.println(stats);
//  }
}
