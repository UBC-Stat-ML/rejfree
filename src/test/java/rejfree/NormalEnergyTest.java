package rejfree;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import bayonet.math.NumericalUtils;



public class NormalEnergyTest
{
  @Test
  public void test()
  {
    double [] means = new double[]{0,0};
    double [][] covar = new double[][]{{1,0.5},{0.5,1}};
    MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means, covar);
    NormalEnergy energy = NormalEnergy.withCovariance(new DoubleMatrix(covar));
    double [] point = new double[]{-0.4, 1.2};
    
    Assert.assertEquals(-energy.valueAt(point),Math.log(mnd.density(point)),NumericalUtils.THRESHOLD);
  }
}


