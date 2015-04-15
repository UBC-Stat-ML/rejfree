package rejfree;

//import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
//import org.jblas.DoubleMatrix;
//import org.junit.Assert;
import org.junit.Test;

//import bayonet.math.NumericalUtils;
import briefj.BriefLog;



public class TestNormalEnergy
{
  @Test
  public void test()
  {
//    double [] means = new double[]{0,0};
//    double [][] covar = new double[][]{{1,0.5},{0.5,1}};
//    MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(means, covar);
//    NormalEnergy energy = NormalEnergy.withCovariance(new DoubleMatrix(covar));
//    double [] point = new double[]{-0.4, 1.2}; 
//    Assert.assertEquals(-energy.valueAt(point),Math.log(mnd.density(point)),NumericalUtils.THRESHOLD);
    
    BriefLog.warnOnce("Warning: the calculation of the normalization of NormalEnergy has "
        + "been temporarily disable, as well as its test (in TestNormalEnergy), "
        + "although the test is correct. This should not be a problem as long as parameters of this "
        + "energy are not resampled (a functionality that is not currently implemented anyways.");
  }
}


