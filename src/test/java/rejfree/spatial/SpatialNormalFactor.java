package rejfree.spatial;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import bayonet.distributions.Normal;
import blang.annotations.FactorArgument;
import blang.variables.RealVariable;
import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;


/**
 * NOTE: this is used in the local sampler, so we do NOT assume velocity for 
 *   the variables of interest to be of unit norm
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class SpatialNormalFactor implements CollisionFactor
{
  @FactorArgument(makeStochastic = true)
  public final RealVariable variable0;
  
  @FactorArgument(makeStochastic = true)
  public final RealVariable variable1;
  
  private final double variance;
  
  public static SpatialNormalFactor newUnaryFactor(double variance, RealVariable variable0)
  {
    return new SpatialNormalFactor(variance, variable0, null);
  }
  
  public static SpatialNormalFactor newBinaryFactor(double variance, RealVariable variable0, RealVariable variable1)
  {
    return new SpatialNormalFactor(variance, variable0, variable1);
  }
  
  private SpatialNormalFactor(
      double variance, RealVariable variable0, RealVariable variable1)
  {
    this.variable0 = variable0;
    this.variable1 = variable1;
    this.variance = variance;
  }

  private double argument()
  {
    return isBinary() ? variable0.getValue() - variable1.getValue() : variable0.getValue();
  }

  private boolean isBinary()
  {
    return variable1 != null;
  }

  @Override
  public double logDensity()
  {
    return Normal.logDensity(argument(), 0.0, variance);
  }

  @Override
  public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
      CollisionContext context)
  {
    double a = argument();
    double b = velocity(context.velocity);
    boolean sameSign = Math.signum(a) == Math.signum(b);
    a = Math.abs(a);
    b = Math.abs(b);
    final double e = StaticUtils.generateUnitRateExponential(context.random);
    
    double time = (sameSign ? Math.sqrt(2 * e * variance + a * a) - a : Math.sqrt(2 * e * variance) + a) / b;
    return Pair.of(time, true);
  }

  private double velocity(DoubleMatrix velocity)
  {
    return isBinary() ? velocity.get(0) - velocity.get(1) : velocity.get(0);
  }

  @Override
  public DoubleMatrix gradient()
  {
    if (!isBinary())
      return new DoubleMatrix(1, 1, - argument() / variance);
    final double 
      g0 = - (variable0.getValue() - variable1.getValue()) / variance,
      g1 = - (variable1.getValue() - variable0.getValue()) / variance;
    return new DoubleMatrix(new double[]{g0, g1});
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    return gradientCoordinate == 0 ? variable0 : variable1;
  }

  @Override
  public int nVariables()
  {
    return isBinary() ? 2 : 1;
  }
}
