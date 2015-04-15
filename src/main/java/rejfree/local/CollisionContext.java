package rejfree.local;

import java.util.Random;

import org.jblas.DoubleMatrix;



public class CollisionContext
{
  public final Random random;
  public final DoubleMatrix velocity;
  public CollisionContext(Random random, DoubleMatrix velocity)
  {
    this.random = random;
    this.velocity = velocity;
  }
}
