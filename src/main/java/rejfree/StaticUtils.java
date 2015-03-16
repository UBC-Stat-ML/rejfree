package rejfree;

import java.util.Random;

import org.jblas.DoubleMatrix;



public class StaticUtils
{
  public static DoubleMatrix position(DoubleMatrix initialPos, DoubleMatrix velocity, double time)
  {
    return initialPos.add(velocity.mul(time));
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
}
