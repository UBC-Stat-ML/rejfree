package rejfree;


import org.jblas.DoubleMatrix;
import org.jblas.util.Random;
import org.junit.Assert;
import org.junit.Test;

import bayonet.math.NumericalUtils;



public class TestBounce
{
  @Test
  public void testProperty1()
  {
    int dim = 3;
    Random.seed(10001);
    
    DoubleMatrix 
      gradient = DoubleMatrix.rand(dim);
    
    DoubleMatrix bounced = StaticUtils.bounce(gradient, gradient);
    Assert.assertEquals(bounced.distance1(gradient.mul(-1.0)), 0.0, NumericalUtils.THRESHOLD);
  }
  
  @Test
  public void testBouncePreservesNorm()
  {
    int dim = 5;
    Random.seed(101);
    
    DoubleMatrix 
      gradient = DoubleMatrix.rand(dim),
      velocityBefore = DoubleMatrix.rand(dim);
    
    double normBefore = velocityBefore.norm2();
    DoubleMatrix bounced = StaticUtils.bounce(velocityBefore, gradient);
    double normAfter = bounced.norm2();
    System.out.println(normBefore);
    System.out.println(normAfter);
    Assert.assertEquals(normBefore, normAfter, NumericalUtils.THRESHOLD);
  }
}
