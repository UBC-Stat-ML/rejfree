package rejfree;

import org.jblas.DoubleMatrix;

import blang.factors.Factor;
import blang.variables.RealVariable;



public interface CollisionFactor extends Factor
{
  public double getCollisionTime(double exponentialRealization, DoubleMatrix velocity);
  public DoubleMatrix gradient();
  public RealVariable getVariable(int gradientCoordinate);
  public int nVariables();
}