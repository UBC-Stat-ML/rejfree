package rejfree.models.normal;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
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
public class RNNormalFactor implements CollisionFactor
{
  @FactorComponent
  public final FactorList<RealVariable> variables;
  
  private final int n;
  private final double delta;
  
  public RNNormalFactor(List<RealVariable> variables)
  {
    this.variables = FactorList.ofArguments(variables, true);
    this.n = variables.size();
    this.delta = 1.0 / ((double) n);
  }
  
  public double getVar(int i)
  {
    double stdDev = ((double) (i+1)) * delta;
    return stdDev * stdDev;
  }

  @Override
  public double logDensity()
  {
    DoubleMatrix point = getPosition();
    return - 0.5 * dotProd(point, point); 
  }
  
  public int dim()
  {
    return n;
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
  
  private double dotProd(DoubleMatrix m1, DoubleMatrix m2)
  {
    double result = 0.0;
    for (int i = 0; i < n; i++)
      result += m1.get(i) * m2.get(i) / getVar(i);
    return result;
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
    
    for (int i = 0; i < n; i++)
      x.put(i, -x.get(i) / getVar(i));
    
    return x;
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    return variables.list.get(gradientCoordinate);
  }

  @Override
  public int nVariables()
  {
    return n;
  }
}