package rejfree;

import java.util.Random;

import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import bayonet.math.NumericalUtils;


public class TestMassMatrix
{
  @Test
  public void testBounceWithMassMatrix()
  {
    org.jblas.util.Random.seed(10001);
    int dim = 10;
    Random rand = new Random(1);
    
    for (int i = 0; i < 10_000; i++)
    {
      DoubleMatrix gradient = DoubleMatrix.randn(dim);
      DoubleMatrix velocity = DoubleMatrix.randn(dim);
      DoubleMatrix massMatrix = randomPosDef(rand, dim);
      
      Assert.assertEquals(
          + dot(gradient, bounce(velocity, gradient, massMatrix), massMatrix),
          - dot(gradient, velocity, massMatrix),
          NumericalUtils.THRESHOLD);
    }
  }
  
  @Test
  public void testNormalTargetInvar()
  {
    org.jblas.util.Random.seed(10001);
    int dim = 10;
    Random rand = new Random(1);
    
    for (int i = 0; i < 10_000; i++)
    {
      DoubleMatrix v = DoubleMatrix.randn(10);
      DoubleMatrix massMatrix = randomPosDef(rand, dim);
      DoubleMatrix gradient = DoubleMatrix.randn(dim);
      
      Assert.assertEquals(
          normalKernel(v, massMatrix),
          normalKernel(bounce(v, gradient, massMatrix), massMatrix),
          NumericalUtils.THRESHOLD);
    }
  }
  
  private static double normalKernel(DoubleMatrix x, DoubleMatrix massMatrix) 
  {
//    DoubleMatrix inv = JBlasUtils.inversePositiveMatrix(massMatrix);
    return dot(x, x, massMatrix);
  }
  
  @Test
  public void testRandomPosDef() 
  {
    org.jblas.util.Random.seed(10001);
    Random rand = new Random(1);
    int dim = 10;
    for (int i = 0; i < 10_000; i++) 
    {
      DoubleMatrix x = DoubleMatrix.randn(dim);
      Assert.assertTrue(dot(x, x, randomPosDef(rand, dim)) > 0.0);
    }
  }
  
  private static DoubleMatrix randomPosDef(Random rand, int dim) 
  {
    DoubleMatrix L = new DoubleMatrix(dim, dim);
    for (int i = 0; i < dim; i++) 
      for (int j = 0; j < dim; j++)
        if (i < j)
          L.put(i, j, rand.nextGaussian());
        else if (i == j) 
          L.put(i, j, rand.nextDouble());
    return L.mmul(L.transpose());
  }
  
  public static DoubleMatrix bounce(DoubleMatrix oldVelocity, DoubleMatrix gradient, DoubleMatrix massMatrix) 
  {
    double scalar = 2.0 * dot(gradient, oldVelocity, massMatrix) / dot(gradient, gradient, massMatrix);
    return oldVelocity.sub(gradient.mul(scalar));
  }
  
  public static double dot(DoubleMatrix v1, DoubleMatrix v2, DoubleMatrix massMatrix) 
  {
    return v1.transpose().mmul(massMatrix).mmul(v2).scalar();
  }
}
