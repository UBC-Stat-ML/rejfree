package rejfree.scalings;

import java.util.List;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import com.google.common.collect.Lists;

public class EffectiveSize2
{
  /**
   * Source: http://www.stats.ox.ac.uk/~burke/Autocorrelation/MCMC%20Output.pdf
   * @param values
   * @return
   */
  public static double effectiveSize(List<Double> values)
  {
    return ((double) values.size())/estimateAutoCorrTime(values);
  }
 
  /**
   * Source: http://www.utstat.toronto.edu/mthompson/MT-ACTMethods.pdf
   * @param values
   * @return
   */
  public static double estimateAutoCorrTime(List<Double> values)
  {
    double sum = 0.0;
    List<Double> autocorrelationFunction = estimateAutocorrelationFunction(values);
    for (int i = 1; i < autocorrelationFunction.size(); i++)
      sum += 2 * autocorrelationFunction.get(i);
    return 1.0 + sum;
  }
  
  public static List<Double> estimateAutocorrelationFunction(List<Double> series)
  {
    SummaryStatistics ss = new SummaryStatistics();
    for (double v : series) ss.addValue(v);
    final double sampleVar = ss.getVariance();
    double [] data = new double[series.size()];
    for (int i =0 ; i < series.size(); i++)
      data[i] = series.get(i);
    double [] auto = autocovariance(data);
    List<Double> result = Lists.newArrayList();
    final double n = series.size();
    double factor = 1.0 / n / sampleVar;
    result.add(1.0);
    loop:for (int i = 1; i < n; i++)
    {
      if (i > 1 && auto[i] + auto[i-1] < 0.0) 
        break loop;
      result.add( factor * auto[i]);
    }
    return result;
  }
  
  /**
   * 
   * Computes the autocovariance of the data in f
   * @param x a vector of real data
   * @param maxShift the maximum phase shift to calculate
   * @return the autocovariance values, having length Math.min(x.length, maxShift)
   * 
   * @author timpalpant
   */
  public static double[] autocovariance(double[] x, int maxShift) {
    double total = 0;
    for (int i = 0; i < x.length; i++) {
      total += x[i];
    }
    double mean = total / x.length;

    int stop = Math.min(x.length, maxShift);
    double[] auto = new double[stop];
    for (int i = 0; i < stop; i++) {
      for (int j = 0; j < x.length - i; j++) {
        auto[i] += (x[j]-mean) * (x[j + i]-mean);
      }
    }

    return auto;
  }

  /**
   * Computes the autocovariance of the data in f for all possible shifts
   * @param x a vector of real data
   * @return the autocovariance values, having length equal to x.length
   * 
   * @author timpalpant
   */
  public static double[] autocovariance(double[] x) {
    return autocovariance(x, x.length);
  }
}
