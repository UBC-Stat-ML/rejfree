package rejfree;

import java.util.List;

import org.jblas.DoubleMatrix;

import blang.variables.RealVariable;



public class NormalFactor implements CollisionFactor
{
  private final List<RealVariable> variables;
  private final NormalEnergy energyFunction;
  private final PegasusConvexCollisionSolver solver;
  
  private NormalFactor(List<RealVariable> variables, NormalEnergy energyFunction)
  {
    this.variables = variables;
    this.energyFunction = energyFunction;
    this.solver = new PegasusConvexCollisionSolver();
  }

  @Override
  public double logDensity()
  {
    return - energyFunction.valueAt(getVector().data);
  }
  
  private DoubleMatrix _vector = null;
  private DoubleMatrix getVector()
  {
    final int dim = variables.size();
    if (_vector == null)
      _vector = new DoubleMatrix(dim);
    for (int i = 0; i < dim; i++)
      _vector.put(i, variables.get(i).getValue());
    return _vector;
  }

  @Override
  public double getCollisionTime(double exponentialRealization, DoubleMatrix velocity)
  {
    // TODO: this could be done analytically
    return solver.collisionTime(getVector(), velocity, energyFunction, exponentialRealization);
  }

  @Override
  public DoubleMatrix gradient()
  {
    DoubleMatrix result = new DoubleMatrix(energyFunction.derivativeAt(getVector().data)).muli(-1.0);
    return result;
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    return variables.get(gradientCoordinate);
  }

  @Override
  public int nVariables()
  {
    return variables.size();
  }

}
