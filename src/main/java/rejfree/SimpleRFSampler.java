package rejfree;

import java.util.Random;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;
import org.jblas.DoubleMatrix;

import bayonet.opt.BacktrackingLineSearcher;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.LBFGSMinimizer;



public class SimpleRFSampler
{
  /**
   * *Negative* log density of the target distribution.
   */
  private final DifferentiableFunction energy;
  
  private final BacktrackingLineSearcher searcher = new BacktrackingLineSearcher();
  private final PegasusSolver solver = new PegasusSolver();
  
  private DoubleMatrix currentPosition, currentVelocity;
  
  public SimpleRFSampler(DifferentiableFunction energy)
  {
    this.energy = energy;
    this.currentPosition = initialPosition(energy);
    this.currentVelocity = initialVelocity(energy.dimension());
  }
  
  private DoubleMatrix initialVelocity(int dimension)
  {
    DoubleMatrix random = DoubleMatrix.randn(dimension);
    double norm = random.norm2();
    return random.muli(1.0/norm);
  }

  private static DoubleMatrix initialPosition(DifferentiableFunction energy)
  {
    double [] min = new LBFGSMinimizer().minimize(energy, new DoubleMatrix(energy.dimension()).data, 1e-2);
    return new DoubleMatrix(min);
  }

  public void iterate(Random rand, int numberOfIterations)
  {
    for (int iter = 0; iter < numberOfIterations; iter++)
    {
      final double collisionTime = collisionTime(currentPosition, currentVelocity, rand.nextDouble());
      currentPosition = position(currentPosition, currentVelocity, collisionTime);
      currentVelocity = Bouncer.bounce(currentVelocity, gradient(currentPosition));
    }
  }

  private DoubleMatrix gradient(DoubleMatrix position)
  {
    return new DoubleMatrix(energy.derivativeAt(position.data));
  }

  private double collisionTime(final DoubleMatrix initialPoint, final DoubleMatrix velocity, final double uniform)
  {
    // go to minimum energy for free
    final DoubleMatrix directionalMin = new DoubleMatrix(searcher.minimize(energy, initialPoint.data, velocity.data));
    final double time1 = time(initialPoint, directionalMin, velocity);
    
    // keep moving until an exponentially distributed amount of energy is exhausted
    final double exponential = - Math.log(uniform);
    final double initialEnergy = energy.valueAt(directionalMin.data);
    final UnivariateFunction lineSolvingFunction = new UnivariateFunction() {
      @Override
      public double value(final double time)
      {
        final DoubleMatrix candidatePosition = position(directionalMin, velocity, time);
        final double candidateEnergy = energy.valueAt(candidatePosition.data);
        final double delta = candidateEnergy - initialEnergy;
        if (delta < 0.0)
          throw new RuntimeException("Did not expect negative delta for convex objective. " +
          		"Delta=" + delta + ", time=" + time);
        return exponential - delta;
      }
    };
    final double upperBound = findUpperBound(lineSolvingFunction);
    final int maxEval = 100;
    final double time2 = solver.solve(maxEval, lineSolvingFunction, 0.0, upperBound);
    return time1 + time2;
  }
  
  private double time(DoubleMatrix initialPos, DoubleMatrix finalPosition, DoubleMatrix velocity)
  {
    final double 
      xInit = initialPos.get(0),
      xFinal= finalPosition.get(0),
      v = velocity.get(0);
    return (xFinal - xInit) / v;
  }

  private DoubleMatrix position(DoubleMatrix initialPos, DoubleMatrix velocity, double time)
  {
    return initialPos.add(velocity.mul(time));
  }
  
  private double findUpperBound(UnivariateFunction lineSolvingFunction)
  {
    double result = 1.0;
    final int maxNIterations = Double.MAX_EXPONENT - 1;
    for (int i = 0; i < maxNIterations; i++)
    {
      if (lineSolvingFunction.value(result) < 0.0)
        return result;
      else
        result *= 2.0;
    }
    throw new RuntimeException();
  }
}
