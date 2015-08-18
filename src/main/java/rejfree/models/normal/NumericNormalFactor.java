package rejfree.models.normal;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import rejfree.StaticUtils;
import rejfree.global.PegasusConvexCollisionSolver;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.variables.RealVariable;
import briefj.BriefLog;


/**
 * An implementation of the normal factor where the collision time is 
 * computed numerically. Mostly for debugging and benchmarking purpose, 
 * in most cases, use NormalFactor instead, which compute the collision 
 * time analytically.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class NumericNormalFactor implements CollisionFactor
{
  @FactorComponent
  public final FactorList<RealVariable> variables;
  
  private final NormalEnergy energyFunction;
  private final PegasusConvexCollisionSolver solver;
  
  public static NumericNormalFactor withPrecision(List<RealVariable> variables, DoubleMatrix precision)
  {
    return new NumericNormalFactor(variables, NormalEnergy.withPrecision(precision));
  }
  
  public static NumericNormalFactor withCovariance(List<RealVariable> variables, DoubleMatrix covar)
  {
    return new NumericNormalFactor(variables, NormalEnergy.withCovariance(covar));
  }
  
  private NumericNormalFactor(List<RealVariable> variables, NormalEnergy energyFunction)
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
