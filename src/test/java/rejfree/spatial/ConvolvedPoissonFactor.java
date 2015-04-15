package rejfree.spatial;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import bayonet.distributions.Bernoulli;
import bayonet.distributions.Exponential;
import bayonet.distributions.Poisson;
import bayonet.opt.DifferentiableFunction;
import blang.annotations.FactorArgument;
import blang.variables.IntegerVariable;
import blang.variables.RealVariable;
import rejfree.PegasusConvexCollisionSolver;
import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;



public class ConvolvedPoissonFactor implements CollisionFactor
{
  @FactorArgument
  public final RealVariable theta0;
  
  @FactorArgument
  public final RealVariable theta1;
  
  @FactorArgument(makeStochastic = true)
  public final IntegerVariable observation;
  
  public ConvolvedPoissonFactor unary(RealVariable theta0, IntegerVariable observation)
  {
    return new ConvolvedPoissonFactor(theta0, null, observation);
  }
  
  public ConvolvedPoissonFactor binary(RealVariable theta0, RealVariable theta1, IntegerVariable observation)
  {
    return new ConvolvedPoissonFactor(theta0, theta1, observation);
  }
  
  public ConvolvedPoissonFactor(RealVariable theta0, RealVariable theta1,
      IntegerVariable observation)
  {
    this.theta0 = theta0;
    this.theta1 = theta1;
    this.observation = observation;
    this.solver = new PegasusConvexCollisionSolver();
  }

  private final PegasusConvexCollisionSolver solver;
  
  private final DifferentiableFunction energyFunction = new DifferentiableFunction() {
    
    @Override
    public double valueAt(double[] thetas)
    {
      return - Poisson.logDensity(observation.getIntegerValue(), lambda(thetas));
    }
    
    private double lambda(double [] thetas)
    {
      return Math.exp(thetas[0]) + (isBinary() ? Math.exp(thetas[1]) : 0);
    }
    
    @Override
    public int dimension()
    {
      return nVariables();
    }
    
    @Override
    public double[] derivativeAt(double[] thetas)
    {
      if (isBinary())
      {
        double [] result = new double[2];
        double factor = ((double)observation.getValue()) / lambda(thetas) - 1.0;
        result[0] = - Math.exp(thetas[0]) * factor;
        result[1] = - Math.exp(thetas[1]) * factor;
        return result;
      }
      else
      {
        double [] result = new double[1];
        result[0] = - ( observation.getValue() - Math.exp(thetas[0]) );
        return result;
      }
    }
  };
  
  private boolean isBinary()
  {
    return theta1 != null;
  }
  
  private double lambda()
  {
    return Math.exp(theta0.getValue()) + (isBinary() ? Math.exp(theta1.getValue()) : 0);
  }

  @Override
  public double logDensity()
  {
    return Poisson.logDensity(observation.getIntegerValue(), lambda());
  }

  
  public Pair<Double,Boolean> getLowerBoundForCollisionDeltaTime(CollisionContext context) 
  {
    if (isBinary())
    {
      // Use the thinning method:
      // compute an upper bound on [0, T]
      final double T = 1.0 + context.random.nextDouble();
      final double 
        v0 = context.velocity.get(0),
        v1 = context.velocity.get(1),
        x0 = theta0.getValue(),
        x1 = theta1.getValue();
      double upperBound = Math.max( // since exp(theta1(t)) + exp(theta2(t)) is convex we can look at end points for maxima
        Math.exp(x0) + Math.exp(x1),
        v0 * Math.exp(x0 + v0 * T) + v1 * Math.exp(x1 + v1 * T));
      double firstEvent = Exponential.generate(context.random, upperBound);
      if (firstEvent > T)
        return Pair.of(T, false);
      double pr =  
          (v0 * Math.exp(x0 + v0 * firstEvent) + v1 * Math.exp(x1 + v1 * firstEvent)) * 
          Math.max(0, 1 - observation.getValue() / (Math.exp(x0 + v0 * firstEvent) + Math.exp(x1 + v1 * firstEvent))/upperBound);
      return Pair.of(firstEvent, Bernoulli.generate(context.random, pr));
        
    }
    else
      return Pair.of(
          solver.collisionTime(getVector(), context.velocity, energyFunction, StaticUtils.generateUnitRateExponential(context.random)),
          true);
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
    switch (gradientCoordinate)
    {
      case 0:  return theta0;
      case 1:  return theta1;
      default: throw new RuntimeException();
    }
  }
  
  private DoubleMatrix _vector = null;
  private DoubleMatrix getVector()
  {
    final int dim = nVariables();
    if (_vector == null)
      _vector = new DoubleMatrix(dim);
    _vector.put(0, theta0.getValue());
    if (isBinary())
      _vector.put(1, theta1.getValue());
    return _vector;
  }

  @Override
  public int nVariables()
  {
    return isBinary() ? 2 : 1;
  }

}
