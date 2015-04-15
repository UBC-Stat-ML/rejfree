package rejfree.local;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import blang.factors.Factor;
import blang.variables.RealVariable;



public interface CollisionFactor extends Factor
{
  /**
   * Computer a lower bound for the next collision time.
   * 
   * @param context The information relevant to making this calculation
   * @return A pair where the first item is a time, which could be either
   *    the collision time (in which case the second item should be true)
   *    or a strict lower bound for the collision time (in which case the
   *    second item should be false).
   */
  public Pair<Double,Boolean> getLowerBoundForCollisionDeltaTime(CollisionContext context);
  public DoubleMatrix gradient();
  public RealVariable getVariable(int gradientCoordinate);
  public int nVariables();
}