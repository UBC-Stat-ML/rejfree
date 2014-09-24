package rejfree;

import org.jblas.DoubleMatrix;



public class Bouncer
{

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
