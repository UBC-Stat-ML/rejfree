package rejfree.mixture;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;
import org.junit.Test;

import bayonet.distributions.Bernoulli;
import bayonet.distributions.Normal;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.GradientValidator;
import blang.annotations.FactorArgument;
import blang.factors.GenerativeFactor;
import blang.variables.RealVariable;
import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;


/**
 * Warning: this was not implemented with the goal of efficiency or
 * generality in mind, but rather to quickly create intuitive plots
 * explaining the execution of the algorithm.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class MixtureFactor implements CollisionFactor, GenerativeFactor
{
  @FactorArgument
  public final RealVariable mean0;
  
  @FactorArgument
  public final RealVariable mean1;
  
  @FactorArgument(makeStochastic = true)
  public final RealVariable observation = RealVariable.real();
  
  public final static double mixtureProp = 0.5;  // hardcoded because of gradient() function
  
  public final static double variance = 1.0; // hardcoded because of gradient() function
  
  public MixtureFactor(RealVariable mean0, RealVariable mean1)
  {
    this.mean0 = mean0;
    this.mean1 = mean1;
  }

  @Override
  public double logDensity()
  {
    return Math.log(
          mixtureProp         * Math.exp(Normal.logDensity(observation.getValue(), mean0.getValue(), variance)) 
        + (1.0 - mixtureProp) * Math.exp(Normal.logDensity(observation.getValue(), mean1.getValue(), variance)));
  }

  @Override
  public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
      CollisionContext context)
  {
    DoubleMatrix velocity = context.velocity;
    double T = 10.0 + context.random.nextDouble();
    double bound = Math.max(evaluateBound(velocity,0), evaluateBound(velocity, T));
    double exp = StaticUtils.generateUnitRateExponential(context.random);
    if (exp > T)
      return Pair.of(T, false);
    
    double objective = Math.max(-gradient(exp, velocity).dot(context.velocity), 0);
    if (objective > bound)
      throw new RuntimeException();
    
    double pr = objective / bound;
    
    return Pair.of(exp, Bernoulli.generate(context.random, pr));
  }
  
  private double evaluateBound(DoubleMatrix velocity, double t)
  {
    double max = Double.NEGATIVE_INFINITY;
    double [] means = means();
    for (int i = 0; i < 2; i++)
    {
      double val = velocity.get(i) * (means[i] + velocity.get(i) * t - observation.getValue());
      if (val > max)
        max = val;
    }
    return Math.max(0, 2.0 * max);
  }
  
  private double [] means()
  {
    return new double[]{mean0.getValue(), mean1.getValue()};
  }
  
  private DoubleMatrix gradient(double deltaT, DoubleMatrix velocity)
  {
    double bu0 = mean0.getValue();
    double bu1 = mean1.getValue();
    DoubleMatrix position = StaticUtils.position(new DoubleMatrix(means()), velocity, deltaT);
    mean0.setValue(position.get(0));
    mean1.setValue(position.get(1));
    DoubleMatrix gradient = gradient();
    mean0.setValue(bu0);
    mean1.setValue(bu1);
    return gradient;
  }

  @Override
  public DoubleMatrix gradient()
  {
    DoubleMatrix result = new DoubleMatrix(2);
    double [] means = means();
    double norm = Math.exp(logDensity());
    for (int i = 0; i < 2 ; i++)
    {
      final double delta = means[i] - observation.getValue();
      result.put(i, - (1.0 / 2.0 / Math.sqrt(2*Math.PI)) * Math.exp( - delta * delta / 2.0) * delta / norm);
    }
    return result;
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    return gradientCoordinate == 0 ? mean0 : mean1;
  }

  @Override
  public int nVariables()
  {
    return 2;
  }

  @Override
  public void generate(Random random)
  {
    double value = Normal.generate(random, (Bernoulli.generate(random, mixtureProp) ? mean0 : mean1).getValue(), variance);
    observation.setValue(value);
  }
  
  @Test
  public void testGradient()
  {
    MixtureFactor f = new MixtureFactor(new RealVariable(1), new RealVariable(2));
    DifferentiableFunction fct = new DifferentiableFunction() {
      
      @Override
      public double valueAt(double[] x)
      {
        f.mean0.setValue(x[0]);
        f.mean1.setValue(x[1]);
        return f.logDensity();
      }
      
      @Override
      public int dimension()
      {
        return 2;
      }
      
      @Override
      public double[] derivativeAt(double[] x)
      {
        f.mean0.setValue(x[0]);
        f.mean1.setValue(x[1]);
        return f.gradient().data;
      }
    };
    Random rand = new Random(1);
    for (int i = 0; i < 10000; i++)
    {
      boolean ok = GradientValidator.isRandomAnalyticDirectionalDerivCloseToNumerical(rand, fct);
      if (!ok)
        throw new RuntimeException();
    }
  }

}
