package rejfree.models.normal;


import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
import blang.annotations.FactorArgument;
import blang.variables.RealVariable;


/**
 * NOTE: this is used in the local sampler, so we do NOT assume velocity for 
 *   the variables of interest to be of unit norm
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class UnivariateNormalFactor implements CollisionFactor
{
  @FactorArgument(makeStochastic = true)
  public final RealVariable variable;
  
  private final double variance;
  
  public UnivariateNormalFactor(RealVariable variable, double variance)
  {
    this.variable = variable;
    this.variance = variance;
  }

  @Override
  public double logDensity()
  {
    double x = variable.getValue();
    return - 0.5 * x * x / variance; 
  }
  
  public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
      CollisionContext context)
  {
    final double x = variable.getValue();
    final double v = context.velocity.get(0);
    
    final double xv = x * v / variance;
    final double vv = v * v / variance;
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
    final DoubleMatrix x = new DoubleMatrix(1);
    
    x.put(0, -variable.getValue() / variance);
    
    return x;
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    return variable;
  }

  @Override
  public int nVariables()
  {
    return 1;
  }
}








