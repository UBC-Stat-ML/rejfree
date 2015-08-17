package rejfree;

import static rejfree.StaticUtils.uniformOnUnitBall;

import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFSampler;
import rejfree.local.NumericNormalFactor;
import rejfree.processors.RecordFullTrajectory;
import bayonet.math.NumericalUtils;
import bayonet.opt.DifferentiableFunction;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;


/**
 * Test the local and global samplers using the analytic formulae
 * available for the normal case.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestAnalyticNormal
{
  /**
   * Note: trajectories should not have refreshments
   * Assume using same sequence of random
   * 
   * @param collisionPoints
   */
  private static void checkAgainsAnalytic(List<DoubleMatrix> collisionPoints)
  {
    Random rand = new Random(1);
    
    // exhaust the first random (used to initialize velocity)
    uniformOnUnitBall(energy.dimension(), rand);
    
    for (int i = 0;
        i < collisionPoints.size() - 2; 
        i++)
    {
      // infer v
      DoubleMatrix 
        v      = getV(collisionPoints.get(i), collisionPoints.get(i+1)),
        vPrime = getV(collisionPoints.get(i+1), collisionPoints.get(i+2));
      DoubleMatrix 
        x      = collisionPoints.get(i),
        xPrime = collisionPoints.get(i+1);
      double e = - Math.log(rand.nextDouble());
      
      double 
        xv      = v     .dot(x),
        xvPrime = vPrime.dot(xPrime);
      
      double s = (xPrime.sub(x)).norm2();
      double analytic_s = analytic_s(e, xv);
      
      double analytic_xvPrime = analytic_xvPrime(e, xv);
      
      Assert.assertEquals(analytic_s, s, NumericalUtils.THRESHOLD);
      Assert.assertEquals(analytic_xvPrime, xvPrime, NumericalUtils.THRESHOLD);
    }
  }
  
  @Test
  public void testLocalRFSamplerWithIsotropicNormal()
  {
    Random rand = new Random(1);
    NormalModel modelSpec = new NormalModel();
    modelSpec.v1.setValue(0.1);
    modelSpec.v2.setValue(0.2);
    ProbabilityModel model = new ProbabilityModel(modelSpec);
    LocalRFSampler sampler = new LocalRFSampler(model , options);
    RecordFullTrajectory rayProcessor = sampler.addRecordFullTrajectoryProcessor();
    sampler.iterate(rand, nIters);
    checkAgainsAnalytic(rayProcessor.samples);
  }
  
  @Test
  public void testSimpleRFSamplerWithIsotropicNormal()
  {
    Random rand = new Random(1);
    DoubleMatrix initialPosition = new DoubleMatrix(new double[]{0.1, 0.2});
    GlobalRFSampler sampler = new GlobalRFSampler(energy, initialPosition, options);
    sampler.iterate(rand, nIters);
    checkAgainsAnalytic(sampler.getTrajectory());
  }
  
  private static class NormalModel
  {
    RealVariable 
      v1 = RealVariable.real(),
      v2 = RealVariable.real();
    
    @DefineFactor
    NumericNormalFactor normalFactor = NumericNormalFactor.withCovariance(Lists.newArrayList(v1, v2), covar);
  }
  
  private static int nIters = 1000;
  private static DoubleMatrix covar = new DoubleMatrix(new double[][]{{0.5,0},{0,0.5}});
  private static DifferentiableFunction energy = NormalEnergy.withCovariance(covar);
  private static RFSamplerOptions options = initOptions();
  
  private static RFSamplerOptions initOptions()
  {
    RFSamplerOptions result = new RFSamplerOptions();
    result.collectRate = 0.0;
    result.refreshRate = 0.0;
    return result;
  }

  private static double analytic_s(double e, double xv)
  {
    return -xv + Math.sqrt(e + (xv > 0 ? xv*xv : 0));
  }

  private static DoubleMatrix getV(
      DoubleMatrix x1,
      DoubleMatrix x2)
  {
    DoubleMatrix diff = x2.sub(x1);
    return diff.div(diff.norm2());
  }

  private static double analytic_xvPrime(double e, double xv)
  {
    return - Math.sqrt(e + (xv > 0 ? xv*xv : 0));
  }
}
