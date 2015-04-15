package rejfree;

import java.util.Random;

import org.jblas.DoubleMatrix;



public class StaticUtils
{
  public static DoubleMatrix position(DoubleMatrix initialPos, DoubleMatrix velocity, double time)
  {
    return initialPos.add(velocity.mul(time));
  }
  
  public static double generateUnitRateExponential(Random random)
  {
    return - Math.log(random.nextDouble());
  }
  
  /**
   * 
   * @param dimension
   * @param rand
   * @return A random vector of unit norm.
   */
  public static DoubleMatrix uniformOnUnitBall(int dimension, Random rand)
  {
    DoubleMatrix random = new DoubleMatrix(dimension);
    for (int i = 0; i < dimension; i++)
      random.data[i] = rand.nextGaussian();
    double norm = random.norm2();
    return random.muli(1.0/norm);
  }

  /**
   * 
   * @param oldVelocity Row vector of velocities before collision
   * @param gradient Row vector of the gradient of the log density at collision
   * @return Row vector of updated velocities
   */
  public static DoubleMatrix bounce(DoubleMatrix oldVelocity, DoubleMatrix gradient)
  {
    final double scale = 2.0 * gradient.dot(oldVelocity) / gradient.dot(gradient);
    return oldVelocity.sub(gradient.mul(scale)); 
  }
}
