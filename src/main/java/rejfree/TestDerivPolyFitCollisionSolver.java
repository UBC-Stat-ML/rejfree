package rejfree;

import java.util.Arrays;
import java.util.function.DoubleFunction;

import org.junit.Assert;
import org.junit.Test;

import bayonet.math.NumericalUtils;
import rejfree.DerivativePolyFitCollisionSolver.PolynomialFit;



public class TestDerivPolyFitCollisionSolver
{
  
  @Test
  public void testQuadFit()
  {
    double 
      trueA = 5.1,
      trueB = 2.1,
      trueC = 1.1;
  
    DoubleFunction<Double> quadFct = x -> trueA * x * x + trueB * x + trueC;
    
    PolynomialFit quadraticFitter = PolynomialFit.quadraticFitter(100.0, quadFct);
    for (int i = 0; i < 10; i++)
    {
      System.out.println(quadraticFitter.toString());
      double [] fit = quadraticFitter.getFit();
      System.out.println(Arrays.toString(fit));
      Assert.assertEquals(trueC, fit[0], NumericalUtils.THRESHOLD);
      Assert.assertEquals(trueB, fit[1], NumericalUtils.THRESHOLD);
      Assert.assertEquals(trueA, fit[2], NumericalUtils.THRESHOLD);
      Assert.assertEquals(quadFct.apply(quadraticFitter.getTestX()), quadraticFitter.evaluateFitAtTestPoint(), NumericalUtils.THRESHOLD);
      quadraticFitter.reduceScale();
    }
  }
  
  @Test
  public void testLinFit()
  {
    double 
      trueSlope = 5.1,
      trueIntercept = 2.1;
    DoubleFunction<Double> linFct = x -> trueSlope * x  + trueIntercept;
    
    PolynomialFit quadraticFitter = PolynomialFit.linearFitter(100.0, linFct);
    for (int i = 0; i < 10; i++)
    {
      System.out.println(quadraticFitter.toString());
      double [] fit = quadraticFitter.getFit();
      System.out.println(Arrays.toString(fit));
      Assert.assertEquals(trueIntercept, fit[0], NumericalUtils.THRESHOLD);
      Assert.assertEquals(trueSlope, fit[1], NumericalUtils.THRESHOLD);
      Assert.assertEquals(linFct.apply(quadraticFitter.getTestX()), quadraticFitter.evaluateFitAtTestPoint(), NumericalUtils.THRESHOLD);
      quadraticFitter.reduceScale();
    }
  }
  
}
