package rejfree;

import org.jblas.DoubleMatrix;

import blang.factors.Factor;
import blang.variables.RealVariable;



public interface CollisionFactor extends Factor
{
  public double getCollisionTime(double exponentialRealization);
  public DoubleMatrix gradient();
  public RealVariable getVariable(int gradientCoordinate);
}