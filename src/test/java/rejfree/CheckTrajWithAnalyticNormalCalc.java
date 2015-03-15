package rejfree;

import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;
import org.junit.Assert;

import rejfree.SimpleRFSampler.SimpleRFSamplerOptions;
import bayonet.math.NumericalUtils;
import bayonet.opt.DifferentiableFunction;



public class CheckTrajWithAnalyticNormalCalc
{
  public static void main(String [] args)
  {
    // Note: trajectories should not have refreshments
    // Assume using same sequence of random
    
    List<DoubleMatrix> collisionPoints = null;
    final int nIters = 1000;
    
    // Numerical:
    {
      Random randOriginal = new Random(1);
      DoubleMatrix covar = new DoubleMatrix(new double[][]{{0.5,0},{0,0.5}});
      DifferentiableFunction energy = NormalEnergy.withCovariance(covar);
      SimpleRFSamplerOptions options = new SimpleRFSamplerOptions();
      options.collectRate = 0.0;
      options.refreshRate = 0.0;
//      DoubleMatrix initPos = new DoubleMatrix(new double[]{0.2,0.4});
      SimpleRFSampler sampler //= new SimpleRFSampler(energy, initPos, options);
         = SimpleRFSampler.initializeRFWithLBFGS(energy, options);
      sampler.iterate(randOriginal, nIters);
      collisionPoints = sampler.getTrajectory();
    }
    
    System.out.println("===");
    
    // Analytic:
    {
      Random rand = new Random(1);
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
        System.out.println(e);
        
        double 
          xv      = v     .dot(x),
          xvPrime = vPrime.dot(xPrime);
        
        double s = (xPrime.sub(x)).norm2();
        double analytic_s = analytic_s(e, xv);
        
        double analytic_xvPrime = analytic_xvPrime(e, xv);
        
        System.out.println("xv = " + xv);
        
        System.out.println("s:\texact=" + analytic_s + ", \tnumerical=" + s);
        System.out.println("xvPrim:\texact=" + analytic_xvPrime + ", \tnumerical=" + xvPrime);
        
        Assert.assertEquals(analytic_s, s, NumericalUtils.THRESHOLD);
        Assert.assertEquals(analytic_xvPrime, xvPrime, NumericalUtils.THRESHOLD);
        
        System.out.println("---");
      }
    }
    
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
