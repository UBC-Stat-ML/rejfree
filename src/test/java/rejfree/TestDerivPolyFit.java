package rejfree;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import bayonet.distributions.Uniform;



public class TestDerivPolyFit
{

  // use normal to compare everything
  @Test
  public void testNormalCollision()
  {
    Random rand = new Random(17);
    for (int i = 0; i < 100; i++)
    {
      // normal factor
      double rho = Uniform.generate(rand, -0.99, 0.99);
      double var1 = 0.1 + rand.nextDouble();
      double var2 = 0.1 + rand.nextDouble();
      double[][] data = new double[][]{{var1, var1 * var2 * rho}, {var1 * var2 * rho, var2}};
      DoubleMatrix covar = new DoubleMatrix(data);
      System.out.println(covar);
      NormalEnergy e = NormalEnergy.withCovariance(covar);
      
      // random direction and point
      DoubleMatrix randomPoint = StaticUtils.uniformOnUnitBall(2, rand).mul(rand.nextGaussian());
      DoubleMatrix randomDirection = StaticUtils.uniformOnUnitBall(2, rand).mul(rand.nextDouble());
      
      long sharedSeed = rand.nextLong();
      
      // method 1 
      PegasusConvexCollisionSolver solver1 = new PegasusConvexCollisionSolver();
      Random method1Rand = new Random(sharedSeed);
      double exp = StaticUtils.generateUnitRateExponential(method1Rand);
      double collisionTime1 = solver1.collisionTime(randomPoint, randomDirection, e, exp);
      System.out.println("gold standard = " + collisionTime1);
      
      // method 2
      DerivativePolyFitCollisionSolver solver2 = new DerivativePolyFitCollisionSolver();
      solver2.initialMax_x = collisionTime1 * 2.0;
      solver2.threshold = 1e-10;
      Random method2Rand = new Random(sharedSeed);
      Pair<Double, Boolean> collision2 = solver2.collisionTime(randomPoint, randomDirection, e, method2Rand);
      System.out.println("tested method = " + collision2);
      
      Assert.assertTrue(collision2.getRight());
      Assert.assertEquals(collisionTime1, collision2.getLeft(), 1e-4);
    }
  }
}
