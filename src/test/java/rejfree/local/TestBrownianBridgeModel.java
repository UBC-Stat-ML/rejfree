package rejfree.local;

import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import rejfree.models.normal.BrownianBridge;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import cern.colt.matrix.tdouble.algo.decomposition.SparseDoubleCholeskyDecomposition;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

import com.google.common.base.Stopwatch;



public class TestBrownianBridgeModel implements Runnable
{
  @OptionSet(name = "runner")
  public LocalRFRunner runner = new LocalRFRunner();
  
  @Option
  public int nVariables = 5;
  
  @Test
  public void testCholesky()
  {
    BrownianBridge model = BrownianBridge.regularlySpaced(nVariables);
    
    SparseDoubleMatrix2D precisionMatrix = model.precisionMatrix();
    SparseDoubleCholeskyDecomposition decomp = new SparseDoubleCholeskyDecomposition(precisionMatrix.getColumnCompressed(false), 0);
    
    DoubleMatrix denseMatrix = new DoubleMatrix(precisionMatrix.toArray());
    DoubleMatrix cholesky = Decompose.cholesky(denseMatrix);
    
    DoubleMatrix sparseL = new DoubleMatrix(decomp.getL().toArray());
    DoubleMatrix denseL = new DoubleMatrix(cholesky.toArray2()).transpose();
    
    System.out.println(sparseL);
    System.out.println(denseL);
    
    for (int i = 0; i < sparseL.rows; i++)
      for (int j = 0; j < sparseL.columns; j++)
        Assert.assertEquals(sparseL.get(i,j), denseL.get(i,j), 0.00001);
  }
  
  @Test
  public void testDeterminant()
  {
    BrownianBridge model = BrownianBridge.regularlySpaced(nVariables);
    
    SparseDoubleMatrix2D precisionMatrix = model.precisionMatrix();
    
    RealMatrix denseMatrix = MatrixUtils.createRealMatrix(precisionMatrix.toArray());
    
    double logAbsDetFromDense  = Math.log(Math.abs(new LUDecomposition(denseMatrix).getDeterminant()));
    double logAbsDetFromSparse = model.logAbsDetPrecision();
    
    System.out.println(logAbsDetFromDense);
    System.out.println(logAbsDetFromSparse);
    
    Assert.assertEquals(logAbsDetFromDense, logAbsDetFromSparse, 1e-10);
  }
  
  @Test
  public void run()
  {
    final BrownianBridge model = BrownianBridge.regularlySpaced(nVariables);
    
    runner.options.maxSteps = 100_000;
    runner.init(model.localFactorModelSpec());
    runner.addMomentRayProcessor();
    runner.run();
    
    for (int i = 0; i < nVariables; i++)
    {
      final double mc = runner.momentRayProcessor.getVarianceEstimate(model.variables.get(i));
      final double analytic = model.getMarginalVariance(i);
          
      System.out.println("MC: " + mc);
      System.out.println("Analytic: " + analytic);
      System.out.println();
      
      Assert.assertEquals(mc, analytic, 0.01);
    }
    
  }
  
  public static void main(String [] args)
  {
    // warmup
    for (int i = 0; i < 10; i++) 
    {
      BrownianBridge bbm = BrownianBridge.regularlySpaced(10);
      bbm.logAbsDetPrecision();
    }
    
    // make sure scales linearly
    for (int nVar = 2; nVar < 1_000_000; nVar *= 2)
    {
      Stopwatch watch = Stopwatch.createStarted();
      BrownianBridge bbm = BrownianBridge.regularlySpaced(nVar);
      bbm.logAbsDetPrecision();
      watch.stop();
      System.out.printf("%d\t%d ms\n", nVar, watch.elapsed(TimeUnit.MILLISECONDS));
    }
    
    Mains.instrumentedRun(args, new TestBrownianBridgeModel());
  }
}
