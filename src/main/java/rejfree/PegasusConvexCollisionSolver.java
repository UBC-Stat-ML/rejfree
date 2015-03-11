package rejfree;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;
import org.jblas.DoubleMatrix;

import bayonet.math.NumericalUtils;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.LBFGSMinimizer;

import static rejfree.StaticUtils.*;

public class PegasusConvexCollisionSolver
{
  private final PegasusSolver solver = new PegasusSolver();
  
  public double collisionTime(final DoubleMatrix initialPoint, final DoubleMatrix velocity, DifferentiableFunction energy, final double exponential)
  {
    // go to minimum energy for free
    final DoubleMatrix directionalMin = lineMinimize(initialPoint, velocity, energy);//new DoubleMatrix(searcher.minimize(energy, initialPoint.data, velocity.data));
    final double time1 = time(initialPoint, directionalMin, velocity);
    
    // keep moving until an exponentially distributed amount of energy is exhausted
    
    final double initialEnergy = energy.valueAt(directionalMin.data);
    final UnivariateFunction lineSolvingFunction = new UnivariateFunction() {
      @Override
      public double value(final double time)
      {
        final DoubleMatrix candidatePosition = position(directionalMin, velocity, time);
        final double candidateEnergy = energy.valueAt(candidatePosition.data);
        final double delta = candidateEnergy - initialEnergy;
        if (delta < - NumericalUtils.THRESHOLD)
          System.err.println("Did not expect negative delta for convex objective. " +
              "Delta=" + delta + ", time=" + time);
        return exponential - delta;
      }
    };
    final double upperBound = findUpperBound(lineSolvingFunction);
    final int maxEval = 100;
    final double time2 = solver.solve(maxEval, lineSolvingFunction, 0.0, upperBound);
    return time1 + time2;
  }
  
  private static DoubleMatrix lineMinimize(
      final DoubleMatrix initialPoint,
      final DoubleMatrix velocity,
      final DifferentiableFunction energy)
  {
    DifferentiableFunction lineRestricted = new DifferentiableFunction() {
      
      @Override
      public double valueAt(double[] _time)
      {
        double time = _time[0];
        double [] position = position(initialPoint, velocity, time).data;
        return energy.valueAt(position);
      }
      
      @Override
      public int dimension()
      {
        return 1;
      }
      
      @Override
      public double[] derivativeAt(double[] _time)
      {
        double time = _time[0];
        double [] position = position(initialPoint, velocity, time).data;
        DoubleMatrix fullDerivative = new DoubleMatrix(energy.derivativeAt(position));
        double directionalDeriv = fullDerivative.dot(velocity);
        return new double[]{directionalDeriv};
      }
    };
    
    double minTime = new LBFGSMinimizer().minimize(lineRestricted, new double[]{0}, 1e-10)[0];
    
    if (minTime < 0.0)
      minTime = 0.0;
    
    return position(initialPoint, velocity, minTime);
  }

  private static double time(DoubleMatrix initialPos, DoubleMatrix finalPosition, DoubleMatrix velocity)
  {
    final double 
      xInit = initialPos.get(0),
      xFinal= finalPosition.get(0),
      v = velocity.get(0);
    return (xFinal - xInit) / v;
  }
  

  
  private static double findUpperBound(UnivariateFunction lineSolvingFunction)
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
