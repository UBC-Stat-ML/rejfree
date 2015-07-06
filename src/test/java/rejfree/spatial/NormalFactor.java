package rejfree.spatial;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;

import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
import bayonet.math.JBlasUtils;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.variables.RealVariable;


/**
 * NOTE: this is used in the local sampler, so we do NOT assume velocity for 
 *   the variables of interest to be of unit norm
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class NormalFactor implements CollisionFactor
{
  @FactorComponent
  public final FactorList<RealVariable> variables;
  
  private final DoubleMatrix precision;
  
  /**
   * log((2pi)^{-k/2} + log(|sigma|^{-1/2}) 
   */
  private final double logConstant; 
  
  public static NormalFactor newUnaryFactor(DoubleMatrix precision, RealVariable variable0)
  {
    List<RealVariable> variables = new ArrayList<RealVariable>();
    variables.add(variable0);
    return new NormalFactor(precision, variables);
  }
  
  public static NormalFactor newBinaryFactor(DoubleMatrix precision, RealVariable variable0, RealVariable variable1)
  {
    List<RealVariable> variables = new ArrayList<RealVariable>();
    variables.add(variable0);
    variables.add(variable1);
    return new NormalFactor(precision, variables);
  }
  
  public static NormalFactor withCovariance(DoubleMatrix covar, List<RealVariable> variables)
  {
    return new NormalFactor(JBlasUtils.inversePositiveMatrix(covar), variables);
  }
  
  private NormalFactor(DoubleMatrix precision, List<RealVariable> variables)
  {
    this.variables = FactorList.ofArguments(variables, true);
    this.precision = precision.dup();
    double detPrecAbs = Math.abs(Decompose.lu(precision).u.diag().prod());
    this.logConstant = - (((double)precision.getRows()) / 2.0) * Math.log(2.0 * Math.PI) + 0.5 * Math.log(detPrecAbs);
  }

  @Override
  public double logDensity()
  {
    DoubleMatrix point = getPosition();
    return - 0.5 * dotProd(point, point) + logConstant; 
  }
  
  private double dotProd(DoubleMatrix x1, DoubleMatrix x2)
  {
    return x1.transpose().mmul(precision).mmul(x2).get(0);
  }
  
  public int dim()
  {
    return precision.rows;
  }

  private DoubleMatrix getPosition()
  {
    DoubleMatrix result = new DoubleMatrix(dim());
    for (int i = 0; i < dim(); i++)
      result.put(i, variables.list.get(i).getValue());
    return result;
  }
  
  public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
      CollisionContext context)
  {
    final DoubleMatrix x = getPosition();
    final DoubleMatrix v = context.velocity;
    
    final double xv = dotProd(x, v);
    final double vv = dotProd(v, v);
    final double e = StaticUtils.generateUnitRateExponential(context.random);
    
    return Pair.of(normalCollisionTime(e, xv, vv),true); 
  }
  
  public static double normalCollisionTime(double exponential, double xv, double vv)
  {
    final double s1 = xv < 0 ? - xv / vv : 0.0;
    final double C = - exponential - s1 * (xv + vv * s1 / 2.0);
    final double result = (- xv + Math.sqrt(xv * xv - 2.0 * vv * C)) / vv;
    return result;
  }

  @Override
  public DoubleMatrix gradient()
  {
    final DoubleMatrix x = getPosition();
    return precision.mmul(x).mul(-1.0);
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    return variables.list.get(gradientCoordinate);
  }

  @Override
  public int nVariables()
  {
    return dim();
  }
}
