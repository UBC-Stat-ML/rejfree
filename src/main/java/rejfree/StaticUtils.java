package rejfree;

import java.util.Random;

import org.jblas.DoubleMatrix;

import bayonet.math.NumericalUtils;



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
  
  /**
   * 
   * @return Array of length 0, 1, or 2 containing the real roots of the polynomial
   *   a x^2 + b x + c
   */
  public static double [] quadraticRealRoots(double a, double b, double c)
  {
    final double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0)
      return new double[0];
    else if (NumericalUtils.isClose(0.0, discriminant, 1e-100))
      return new double[]{ -b / 2.0 / a};
    else
    {
      final double sqrtDiscr = Math.sqrt(discriminant);
      final double twoA = 2.0 * a;
      return new double[]{ 
        (-b - sqrtDiscr) / twoA, 
        (-b + sqrtDiscr) / twoA };
    }
  }
  
  public static double [] quadraticRealRoots(double [] coefficients)
  {
    return quadraticRealRoots(coefficients[2], coefficients[1], coefficients[0]);
  }
  
  public static double indefiniteIntegralOfQuadratic(double [] coefficients, double x)
  {
    final double 
      c = coefficients[0],
      b = coefficients[1],
      a = coefficients[2];
    return a * x * x * x / 3.0 + b * x * x / 2.0 + c * x;
  }
  
  public static double indefiniteIntegralOfLinear(double slope, double intercept, double x)
  {
    return slope * x * x / 2.0 + intercept * x;
  }
  
  public static double indefiniteIntegralOfLinear(double [] coefficients, double x)
  {
    return indefiniteIntegralOfLinear(coefficients[1], coefficients[0], x);
  }
  
  public static double linearRoot(double slope, double intercept)
  {
    return - intercept / slope;
  }
  
  public static double linearRoot(double [] coefficients)
  {
    return linearRoot(coefficients[1], coefficients[0]);
  }
  
  private StaticUtils() {}
}
