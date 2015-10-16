package rejfree.models.expfam;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import bayonet.distributions.Poisson;
import blang.annotations.FactorArgument;
import blang.factors.GenerativeFactor;
import blang.variables.IntegerVariable;
import blang.variables.RealVariable;
import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;



public class PoissonFactor implements CollisionFactor, GenerativeFactor
{
  @FactorArgument
  public final RealVariable naturalParameter;
  
  @FactorArgument(makeStochastic = true)
  public final IntegerVariable realization;
  
  public PoissonFactor(
      RealVariable naturalParameter,
      IntegerVariable realization)
  {
    this.naturalParameter = naturalParameter;
    this.realization = realization;
  }

  private int getRealization()
  {
    return realization.getIntegerValue();
  }
  
  private double mean()
  {
    return Math.exp(naturalParameter.getValue());
  }

  @Override
  public double logDensity()
  {
    return Poisson.logDensity(getRealization(), mean());
  }

  @Override
  public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
      CollisionContext context)
  {
    final double v = context.velocity.get(0);
    final double theta0 = naturalParameter.getValue();
    final int x = getRealization();
    
    double t1 = - StaticUtils.generateUnitRateExponential(context.random) / v / x;
    double t2 = 
      (
        Math.log
        (
          StaticUtils.generateUnitRateExponential(context.random) + Math.exp(theta0)
        ) 
        - theta0
      )
      / v;
    
    t1 = t1 >= 0 ? t1 : Double.POSITIVE_INFINITY;
    t2 = t2 >= 0 ? t2 : Double.POSITIVE_INFINITY;
    
    return Pair.of(Math.min(t1, t2), true);
  }

  @Override
  public DoubleMatrix gradient()
  {
    DoubleMatrix result = new DoubleMatrix(1);
    result.put(0, ((double) getRealization()) - Math.exp(naturalParameter.getValue()));
    return result;
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    if (gradientCoordinate != 0)
      throw new RuntimeException();
    return naturalParameter;
  }

  @Override
  public int nVariables()
  {
    return 1;
  }

  @Override
  public void generate(Random random)
  {
    realization.setValue(Poisson.generate(random, mean()));
  }
}
