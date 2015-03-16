package rejfree;

import static rejfree.StaticUtils.uniformOnUnitBall;

import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;

import rejfree.LocalRFSampler.RecordFullTrajectory;
import rejfree.SimpleRFSampler.RFSamplerOptions;
import bayonet.math.NumericalUtils;
import bayonet.opt.DifferentiableFunction;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;



public class CheckTrajWithAnalyticNormalCalc
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
      
//      System.out.println("xv = " + xv);
      
//      System.out.println("s:\texact=" + analytic_s + ", \tnumerical=" + s);
//      System.out.println("xvPrim:\texact=" + analytic_xvPrime + ", \tnumerical=" + xvPrime);
      
      Assert.assertEquals(analytic_s, s, NumericalUtils.THRESHOLD);
      Assert.assertEquals(analytic_xvPrime, xvPrime, NumericalUtils.THRESHOLD);
      
//      System.out.println("---");
    }
  }
  
  @Test
  public void testLocalRFSamplerWithIsotropicNormal()
  {
    Random rand = new Random(1);
    ProbabilityModel model = new ProbabilityModel(new NormalModel());
    LocalRFSampler sampler = new LocalRFSampler(model , options);
    RecordFullTrajectory rayProcessor = sampler.addRecordFullTrajectoryProcessor();
    sampler.iterate(rand, nIters);
    checkAgainsAnalytic(rayProcessor.samples);
  }
  
  @Test
  public void testSimpleRFSamplerWithIsotropicNormal()
  {
    Random rand = new Random(1);
    SimpleRFSampler sampler = SimpleRFSampler.initializeRFWithLBFGS(energy, options);
    sampler.iterate(rand, nIters);
    checkAgainsAnalytic(sampler.getTrajectory());
  }
  
  private static class NormalModel
  {
    RealVariable 
      v1 = RealVariable.real(),
      v2 = RealVariable.real();
    
    @DefineFactor
    NormalFactor normalFactor = NormalFactor.withCovariance(Lists.newArrayList(v1, v2), covar);
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
