package rejfree.spatial;


import java.util.Random;

import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import bayonet.opt.GradientValidator;
import blang.variables.IntegerVariable;
import blang.variables.RealVariable;
import rejfree.NormalEnergy;



public class TestGradient
{
  @Test
  public void testNormalEnergyGradient()
  {
    DoubleMatrix covar = new DoubleMatrix(new double[][]{{4.1, -0.2},{-0.2, 6.1}});
    NormalEnergy e = NormalEnergy.withCovariance(covar);
    
    Random random = new Random(1);
    for (int i = 0; i < 10; i++)
      Assert.assertTrue(GradientValidator.isRandomAnalyticDirectionalDerivCloseToNumerical(random, e));
  }
  
  @Test
  public void testUnaryPoissonGradient()
  {
    Random rand = new Random(1);
    ConvolvedPoissonFactor cpf = ConvolvedPoissonFactor.unary(new RealVariable(rand.nextDouble()), new IntegerVariable(rand.nextInt(10)));
    
    for (int i = 0; i < 10; i++)
      Assert.assertTrue(GradientValidator.isRandomAnalyticDirectionalDerivCloseToNumerical(rand, cpf.energyFunction));
  }
  
  @Test
  public void testBinaryPoissonGradient()
  {
    Random rand = new Random(1);
    ConvolvedPoissonFactor cpf = ConvolvedPoissonFactor.binary(new RealVariable(rand.nextDouble()), new RealVariable(rand.nextDouble()), new IntegerVariable(rand.nextInt(10)));
    
    for (int i = 0; i < 10; i++)
      Assert.assertTrue(GradientValidator.isRandomAnalyticDirectionalDerivCloseToNumerical(rand, cpf.energyFunction));
  }
  
//  @Test
//  public void testSpatialNormalFactorGradient()
//  {
//    Random rand = new Random(1);
//    double variance = 1.7;
//    RealVariable variable0 = new RealVariable(rand.nextDouble());
//    RealVariable variable1 = new RealVariable(rand.nextDouble());
//    SpatialNormalFactor f = SpatialNormalFactor.newBinaryFactor(variance, variable0, variable1);
//    
//    DifferentiableFunction df = new DifferentiableFunction() {
//      
//      @Override
//      public double valueAt(double[] x)
//      {
//        variable0.setValue(x[0]);
//        variable1.setValue(x[1]);
//        double logDensity = f.logDensity();
//        return logDensity;
//      }
//      
//      @Override
//      public int dimension()
//      {
//        return 2;
//      }
//      
//      @Override
//      public double[] derivativeAt(double[] x)
//      {
//        variable0.setValue(x[0]);
//        variable1.setValue(x[1]);
//        return f.gradient().data;
//      }
//    };
//    
//    for (int i = 0; i < 10; i++)
//      Assert.assertTrue(GradientValidator.isRandomAnalyticDirectionalDerivCloseToNumerical(rand, df));
//  }
  
//  @Test
//  public void testSpatialNormalUnaryFactorGradient()
//  {
//    Random rand = new Random(7);
//    double variance = 1.7;
//    RealVariable variable0 = new RealVariable(rand.nextDouble());
//    SpatialNormalFactor f = SpatialNormalFactor.newUnaryFactor(variance, variable0);
//    
//    DifferentiableFunction df = new DifferentiableFunction() {
//      
//      @Override
//      public double valueAt(double[] x)
//      {
//        variable0.setValue(x[0]);
//        double logDensity = f.logDensity();
//        return logDensity;
//      }
//      
//      @Override
//      public int dimension()
//      {
//        return 1;
//      }
//      
//      @Override
//      public double[] derivativeAt(double[] x)
//      {
//        variable0.setValue(x[0]);
//        return f.gradient().data;
//      }
//    };
//    
//    for (int i = 0; i < 10; i++)
//      Assert.assertTrue(GradientValidator.isRandomAnalyticDirectionalDerivCloseToNumerical(rand, df));
//  }
}
