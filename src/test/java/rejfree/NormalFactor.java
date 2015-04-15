package rejfree;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.variables.RealVariable;
import briefj.BriefLog;



public class NormalFactor implements CollisionFactor
{
  @FactorComponent
  public final FactorList<RealVariable> variables;
  
  private final NormalEnergy energyFunction;
  private final PegasusConvexCollisionSolver solver;
  
  public static NormalFactor withPrecision(List<RealVariable> variables, DoubleMatrix precision)
  {
    return new NormalFactor(variables, NormalEnergy.withPrecision(precision));
  }
  
  public static NormalFactor withCovariance(List<RealVariable> variables, DoubleMatrix covar)
  {
    return new NormalFactor(variables, NormalEnergy.withCovariance(covar));
  }
  
  private NormalFactor(List<RealVariable> variables, NormalEnergy energyFunction)
  {
    BriefLog.warnOnce("WARNING: NormalFactor collision could be computed analytically.");
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
  public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
      CollisionContext context)
  {
    double exponentialRealization = StaticUtils.generateUnitRateExponential(context.random);
    double collisionTime = solver.collisionTime(getVector(), context.velocity, energyFunction, exponentialRealization);
    return Pair.of(collisionTime, true);
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
