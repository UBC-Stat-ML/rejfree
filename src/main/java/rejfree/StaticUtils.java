package rejfree;

import org.jblas.DoubleMatrix;



public class StaticUtils
{
  public static DoubleMatrix position(DoubleMatrix initialPos, DoubleMatrix velocity, double time)
  {
    return initialPos.add(velocity.mul(time));
  }
}
