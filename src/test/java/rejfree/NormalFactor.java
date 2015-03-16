package rejfree;

import java.util.List;

import org.jblas.DoubleMatrix;

import rejfree.local.CollisionFactor;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.variables.RealVariable;



public class NormalFactor implements CollisionFactor
{
  @FactorComponent
  public final FactorList<RealVariable> variables;
  
  private final NormalEnergy energyFunction;
  private final PegasusConvexCollisionSolver solver;
  
  public static NormalFactor withCovariance(List<RealVariable> variables, DoubleMatrix covar)
  {
    return new NormalFactor(variables, NormalEnergy.withCovariance(covar));
  }
  
  private NormalFactor(List<RealVariable> variables, NormalEnergy energyFunction)
  {
    this.variables = FactorList.ofArguments(variables, true);
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
    final int dim = variables.list.size();
    if (_vector == null)
      _vector = new DoubleMatrix(dim);
    for (int i = 0; i < dim; i++)
      _vector.put(i, variables.list.get(i).getValue());
    return _vector;
  }

  @Override
  public double getCollisionDeltaTime(double exponentialRealization, DoubleMatrix velocity)
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
    return variables.list.get(gradientCoordinate);
  }

  @Override
  public int nVariables()
  {
    return variables.list.size();
  }

}
