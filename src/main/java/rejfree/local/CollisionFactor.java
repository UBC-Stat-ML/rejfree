package rejfree.local;

import org.jblas.DoubleMatrix;

import blang.factors.Factor;
import blang.variables.RealVariable;



public interface CollisionFactor extends Factor
{
  public double getCollisionDeltaTime(double exponentialRealization, DoubleMatrix velocity);
  public DoubleMatrix gradient();
  public RealVariable getVariable(int gradientCoordinate);
  public int nVariables();
}