package rejfree;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import bayonet.distributions.Poisson;
import bayonet.distributions.Uniform;

import com.google.common.collect.Lists;



public class Example1D
{
  private static int nIters = 100;
  
  
  public static List<Double> method1(Random rand)
  {
    List<Double> result = Lists.newArrayList();
    for (int i = 0; i < nIters; i++)
    {
      double u = rand.nextDouble();
      double t = Math.sqrt(-2.0 * Math.log(u));
      int pp = Poisson.generate(rand, t);
      for (int j = 0; j < pp; j++)
        result.add(Uniform.generate(rand, -t, t));
    }
    
    return result;
  }
  
  public static List<Double> method2(Random rand)
  {
    List<Double> result = Lists.newArrayList();
    
    for (int i = 0; i < nIters; i++)
      result.add((rand.nextGaussian()));
    
    return result;
  }
  

  public static void main(String [] args)
  {
    Random rand = new Random(1);
    summarize(method1(rand));
    summarize(method2(rand));
  }
  
  public static void summarize(List<Double> samples)
  {
    SummaryStatistics stats = new SummaryStatistics();
    for (double d : samples)
      stats.addValue(d);
    System.out.println(stats);
  }
}
