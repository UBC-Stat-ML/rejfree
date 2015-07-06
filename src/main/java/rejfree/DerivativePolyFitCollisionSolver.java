package rejfree;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;
import java.util.function.DoubleFunction;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import bayonet.distributions.Bernoulli;
import bayonet.math.NumericalUtils;
import bayonet.opt.DifferentiableFunction;
import static rejfree.StaticUtils.*;


/**
 * Uses a polynomial fit to the directional derivative in order to approximate
 * a lower bound for collision. 
 * 
 * Warning: not extensively tested, potential bug
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class DerivativePolyFitCollisionSolver
{
  private final boolean useQuadratic = false;
  double initialMax_x = 1.0;
  private SummaryStatistics max_xSummary = new SummaryStatistics();
  double threshold = 0.1;
  
  public Pair<Double,Boolean> collisionTime(
      final DoubleMatrix initialPoint, 
      final DoubleMatrix velocity, 
      DifferentiableFunction energy, 
      Random random) 
  {
    // 1- build a new function using directional derivatives
    DoubleFunction<Double> directionalDerivative = directionalDerivative(initialPoint, velocity, energy);
    
    // 2- create the fit
    PolynomialFit polyFit = useQuadratic ? 
      PolynomialFit.quadraticFitter(initialMax_x, directionalDerivative) :
      PolynomialFit.linearFitter(initialMax_x, directionalDerivative);  
    
    // 3- rescale until fit is good enough
    while (!satisfactory(polyFit))
      polyFit.reduceScale();
    
    // 4- update stats to make initial scale better next time
    updateScaleStats(polyFit);
    
    // 5a- add fudge factor
    addFudgeFactor(polyFit);
    
    // 5b- integrate
    double bound = integratePositiveLinearPart(polyFit.fit[1], polyFit.fit[0], polyFit.max_x);
    
    // 6- find collision
    double exponential = StaticUtils.generateUnitRateExponential(random);
    
    if (bound < exponential) // no collision in this case:
      return Pair.of(polyFit.max_x, false);
    
    double proposedCollisionTime = solveCollisionTime(polyFit, exponential);
    
    {
      double check = integratePositiveLinearPart(polyFit.fit[1], polyFit.fit[0], proposedCollisionTime);
      System.out.println("" + exponential + " vs " + check);
    }
    
    
    double valueOfBound = polyFit.evaluateFitAt(proposedCollisionTime);
    double valueOfFct = directionalDerivative.apply(proposedCollisionTime);
    
    if (valueOfBound < valueOfFct)
    {
      updateThreshold();
      return collisionTime(initialPoint, velocity, energy, random);
    }
    
    double acceptPr = valueOfBound / valueOfFct;
    boolean accepted = Bernoulli.generate(random, acceptPr);
    
    return Pair.of(proposedCollisionTime, accepted);
  }
  
  private void updateThreshold()
  {
    threshold = 0.9 * threshold;
    max_xSummary = new SummaryStatistics();
    System.err.println("Threshold updated to " + threshold);
  }

  private void updateScaleStats(PolynomialFit polyFit)
  {
    max_xSummary.addValue(polyFit.max_x);
    initialMax_x = 2.0 * max_xSummary.getMean();
  }

  private void addFudgeFactor(PolynomialFit polyFit)
  {
    polyFit.fit[0] += threshold * 2.0;
  }

  private boolean satisfactory(PolynomialFit polyFit)
  {
    double pred  = polyFit.evaluateFitAtTestPoint();
    double truth = polyFit.evaluateFunctionAtTestPoint();
    return NumericalUtils.isClose(pred, truth, threshold);
  }

  private double solveCollisionTime(PolynomialFit polyFit, double exponential)
  {
    if (polyFit.order == 2)
      throw new RuntimeException();
    
    final double
      b = polyFit.fit[0],
      a = polyFit.fit[1];
    
    if (NumericalUtils.isClose(a, 0.0, 10-10))
      return exponential / b;
    else
    {
      final double r = linearRoot(a, b);
      // in the following, we are using the fact that we know a collision will occur
      final double c = - exponential - indefiniteIntegralOfLinear(a, b, a > 0 ? Math.max(0,r) : 0); 
      double[] quadraticRealRoots = quadraticRealRoots(a / 2.0, b, c);
      if (quadraticRealRoots[0] >= r)
        return quadraticRealRoots[0];
      else
        return quadraticRealRoots[1];
      
//      {
//        System.out.println(exponential);
//        System.out.println("--->" + indefiniteIntegralOfLinear(a, b, a > 0 ? Math.max(0,r) : 0));
//        PlotLine.fromFunction(x -> a * x * x / 2.0 + b * x - indefiniteIntegralOfLinear(a, b, a > 0 ? Math.max(0,r) : 0), 0.0, polyFit.max_x).openTemporaryPDF();
//      }
      
    }
  }

//  private double smallestPositiveSolution(double[] quadraticRealRoots)
//  {
//    for (double d : quadraticRealRoots)
//      if (d >= 0)
//        return d;
//    throw new RuntimeException();
//  }

  private DoubleFunction<Double> directionalDerivative(
      final DoubleMatrix initialPoint, 
      final DoubleMatrix velocity,
      final DifferentiableFunction energy)
  {
    return t -> {
      DoubleMatrix point = velocity.mul(t).addi(initialPoint);
      DoubleMatrix gradient = new DoubleMatrix(energy.derivativeAt(point.data));
      return gradient.dot(velocity);
    };
  }

//  private static double integratePositivePart(PolynomialFit polyFit)
//  {
//    if (polyFit.order == 2)
//      return _integratePositiveQuadraticPart(polyFit);
//    else if (polyFit.order == 1)  
//      return _integratePositiveLinearPart(polyFit);
//    else
//      throw new RuntimeException();
//  }
  
//  private static double _integratePositiveQuadraticPart(double [] fit, double T)
//  {
//    final double a = fit[2];
//    
//    if (NumericalUtils.isClose(a, 0.0, 10-10))
//      return _integratePositiveLinearPart(fit, T);
//    else
//    {
//      final double [] roots = quadraticRealRoots(fit);
//      if (roots.length == 2)
//      {
//        final double 
//          r0 = roots[0],
//          r1 = roots[1];
//        if (a > 0) 
//          return indefiniteIntegralOfQuadratic(fit, Math.max(r0, 0)) - indefiniteIntegralOfQuadratic(fit, 0) 
//               + indefiniteIntegralOfQuadratic(fit, T) - indefiniteIntegralOfQuadratic(fit, Math.min(r1, T));
//        else       
//          return indefiniteIntegralOfQuadratic(fit, Math.min(r1, T)) - indefiniteIntegralOfQuadratic(fit, Math.max(r0, 0));
//      }
//      else // no root or single root
//      {
//        if (a > 0) return indefiniteIntegralOfQuadratic(fit, T) - indefiniteIntegralOfQuadratic(fit, 0);
//        else       return 0.0;
//      }
//    }
//  }
  
  private static double integratePositiveLinearPart(double slope, double intercept, double T)
  {
    if (NumericalUtils.isClose(slope, 0.0, 1e-10))
      return T * Math.min(0, intercept);
    else 
    {
      double r = linearRoot(slope, intercept);
      if (slope > 0)  return indefiniteIntegralOfLinear(slope, intercept, Math.max(T, r)) - indefiniteIntegralOfLinear(slope, intercept, Math.max(0, r));
      else            return indefiniteIntegralOfLinear(slope, intercept, Math.min(T, r)) - indefiniteIntegralOfLinear(slope, intercept, Math.min(0, r));
    }
  }
  
  static class PolynomialFit
  {
    private final int order; // 1 for linear, 2 for quadratic
    private final double evaluationAtZero;
    // assumed to be size order + 1 (order+1 to fit, 1 to test) - 1 because of point at zero
    private final LinkedList<Double> additionalEvaluations;
    private double max_x;
    private final double xDecreaseFactor;
    private final DoubleFunction<Double> function;
    private double [] fit;
    
    static PolynomialFit quadraticFitter(
        double max_x, 
        DoubleFunction<Double> function)
    {
      return new PolynomialFit(max_x, function, 2.0, 2);
    }
    
    static PolynomialFit linearFitter(
        double max_x, 
        DoubleFunction<Double> function)
    {
      return new PolynomialFit(max_x, function, 2.0, 1);
    }
    
    private PolynomialFit(
        double max_x, 
        DoubleFunction<Double> function,
        double xDecreaseFactor,
        int order)
    {
      this.order = order;
      this.max_x = max_x;
      this.xDecreaseFactor = xDecreaseFactor;
      this.function = function;
      this.additionalEvaluations = new LinkedList<>();
      double evalAtZero = Double.NaN;
      for (int i = 0; i < order + 2; i++)
      {
        double x = getX(i);
        double y = function.apply(x);
        if (i == 0)
          evalAtZero = y;
        else
          this.additionalEvaluations.add(y);
      }
      this.evaluationAtZero = evalAtZero;
      computeFit();
    }
    
    double getX(int index)
    {
      if (index == 0)
        return 0.0;
      return max_x * Math.pow(xDecreaseFactor, -(index-1));
    }
    
    int nPoints()
    {
      return additionalEvaluations.size() + 1;
    }
    
    void reduceScale()
    {
      // create new x's
      max_x = max_x / xDecreaseFactor;
      // then:
      {
        // shift evals
        additionalEvaluations.remove(0);
        additionalEvaluations.add(function.apply(getTestX()));
      }
      computeFit();
    }
    
    double [] getFit()
    {
      return fit;
    }
    
    double getTestX()
    {
      return getX(order + 1);
    }
    
    double evaluateFunctionAtTestPoint()
    {
      return getY(order + 1);
    }
    
    double getY(int index)
    {
      if (index == 0)
        return evaluationAtZero;
      return additionalEvaluations.get(index - 1);
    }
    
    // TODO: reuse arrays, etc
    private void computeFit()
    {
      if (order == 1)
      {
        double b = getY(0);
        double a = (getY(1) - b)/getX(1);
        fit = new double[]{b, a};
      }
      else if (order == 2)
      {
        // TODO: speedup with a direct solution
        DoubleMatrix coefficients = new DoubleMatrix(new double[][]{
          {1.0, getX(0), getX(0) * getX(0)},
          {1.0, getX(1), getX(1) * getX(1)},
          {1.0, getX(2), getX(2) * getX(2)}
        });
        DoubleMatrix intercepts = new DoubleMatrix(new double[]{
          getY(0),
          getY(1),
          getY(2)
        });
        fit = Solve.solve(coefficients, intercepts).data;
      }
      else throw new RuntimeException();
    }
    
    double evaluateFitAt(double x)
    {
      double result = fit[0] + fit[1] * x;
      if (order == 2)
        result += fit[2] * x * x;
      return result;
    }
    
    double evaluateFitAtTestPoint()
    {
      double x = getTestX();
      return evaluateFitAt(x);
    }
    
    @Override
    public String toString()
    {
      String result = 
        String.format("order = %d%n" + 
          "fit = %s%n", order, Arrays.toString(fit));
      for (int i = 0; i < nPoints(); i++)
        result += String.format("\tx = %f, y = %f%n", getX(i), getY(i));
      return result;
    }
  }
}
