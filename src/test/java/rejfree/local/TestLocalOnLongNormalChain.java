package rejfree.local;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.CollisionFactor;
import rejfree.local.LocalRFSampler;
import bayonet.math.JBlasUtils;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.RealVariable;



public class TestLocalOnLongNormalChain implements Runnable
{
  double diag = 4.0/3.0;
  double offDiag = - 2.0/3.0;
  int nPairs = 4;
  Random random = new Random(1);
  
  private int dim() { return nPairs + 1; }
  
  DoubleMatrix fullPrecision, covarMatrix;
  List<DoubleMatrix> pairPrecisions;
  MultivariateNormalDistribution normal;
  
  private DoubleMatrix exactSample()
  {
    return new DoubleMatrix(normal.sample());
  }
  
  public static double maxError(DoubleMatrix guess, DoubleMatrix truth)
  {
    return (guess.sub(truth)).normmax();
  }
  
  
  private void buildPrecisionMatrices()
  {
    fullPrecision = new DoubleMatrix(nPairs+1, nPairs+1);
    pairPrecisions = new ArrayList<DoubleMatrix>();
    for (int i = 0; i < nPairs; i++)
    {
      DoubleMatrix cur = new DoubleMatrix(2,2);
      pairPrecisions.add(cur);
      for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++)
        {
          double val = r == c ? diag : offDiag;
          cur.put(r, c, val);
          fullPrecision.put(r+i,c+i,val + fullPrecision.get(r+i,c+i));
        }
    }
    covarMatrix = JBlasUtils.inversePositiveMatrix(fullPrecision);
    
    System.out.println("Full precision");
    System.out.println(fullPrecision);
    System.out.println("Full covariance");
    System.out.println(JBlasUtils.inversePositiveMatrix(fullPrecision));
    System.out.println("Pairwise precisions");
    System.out.println(pairPrecisions);
    
    normal = new MultivariateNormalDistribution(new DoubleMatrix(dim()).data, JBlasUtils.asDoubleArray(covarMatrix));
    normal.reseedRandomGenerator(random.nextLong());
  }
  
  public void testMCAvgs(boolean useAnalytic)
  {
    RFSamplerOptions options = new RFSamplerOptions();
    options.refreshRate = 0.0;
    options.collectRate = 1.0;
    final DoubleMatrix 
      outerSumsRF = new DoubleMatrix(dim(),dim());
    int [] i = new int[1];
    DoubleMatrix exactSample = exactSample();
    Model modelSpec = new Model(exactSample.data, useAnalytic);
    ProbabilityModel model = new ProbabilityModel(modelSpec);
    System.out.println(model);
    LocalRFSampler local = new LocalRFSampler(model, options);
    local.addPointProcessor(new Processor()
    {
      @Override
      public void process(ProcessorContext context)
      {
        DoubleMatrix x = currentPositionToVector(modelSpec);
        outerSumsRF.addi(x.mmul(x.transpose()));
        i[0]++;
      }
    });
    local.iterate(random, 100000, Double.POSITIVE_INFINITY);
    outerSumsRF.divi(i[0]);
    System.out.println("empirical from MC avg (nSamples = " + i[0] + ")");
    System.out.println(outerSumsRF);
    
    double maxErr = maxError(covarMatrix, outerSumsRF);
    System.out.println("maxError = " + maxErr);
    Assert.assertTrue(maxErr < 0.01);
  }
  
  public void testInvariance(boolean useAnalytic)
  {
    RFSamplerOptions options = new RFSamplerOptions();
    options.refreshRate = 0.0;
    options.collectRate = 0.0;

    double fixedTime = 2;
    int nRepeats = 20000;
    
    DoubleMatrix 
      outerSumsExact = new DoubleMatrix(dim(),dim()),
      outerSumsRF = new DoubleMatrix(dim(),dim());
    
    for (int i = 0; i < nRepeats; i++)
    {
      DoubleMatrix exactSample = exactSample();
      outerSumsExact.addi(exactSample.mmul(exactSample.transpose()));
      
      DoubleMatrix rfSampl = null;
      Model modelSpec = new Model(exactSample.data, useAnalytic);
      ProbabilityModel model = new ProbabilityModel(modelSpec);
      LocalRFSampler local = new LocalRFSampler(model, options);
      local.iterate(random, Integer.MAX_VALUE, fixedTime);
      rfSampl = currentPositionToVector(modelSpec);
      outerSumsRF.addi(rfSampl.mmul(rfSampl.transpose()));
    }
    outerSumsExact.divi(nRepeats);
    outerSumsRF.divi(nRepeats);
    
    System.out.println("Empirical covar from exact samples");
    System.out.println(outerSumsExact);
    
    System.out.println("Empirical covar from RF samples");
    System.out.println(outerSumsRF);
    
    double maxErr = maxError(covarMatrix, outerSumsRF);
    System.out.println("maxError = " + maxErr);
    Assert.assertTrue(maxErr < 0.02);
  }
  
  static DoubleMatrix currentPositionToVector(Model model)
  {
    int nVars = model.variables.size();
    DoubleMatrix current = new DoubleMatrix(nVars);
    for (int i = 0; i < nVars; i++)
      current.put(i, model.variables.get(i).getValue());
    return current;
  }
  
  @Test
  public void run()
  {
    org.jblas.util.Random.seed(random.nextLong());
    
    buildPrecisionMatrices();
    testMCAvgs(true);
    testInvariance(true);
  }

  public static void main(String[] args)
  {
    new TestLocalOnLongNormalChain().run();
  }
  
  public class Model
  {
    List<RealVariable> variables = new ArrayList<>();
    
    @DefineFactor
    public final List<CollisionFactor> localFactors;
    
    public Model(double [] init, boolean useAnalytic)
    {
      this.localFactors = localFactors(init, useAnalytic);
    }

    private List<CollisionFactor> localFactors(double [] init, boolean useAnalytic)
    {
      List<CollisionFactor> result = new ArrayList<>();
      RealVariable prev = RealVariable.real(init[0]);
      variables.add(prev);
      for (int i = 0; i < init.length - 1; i++)
      {
        DoubleMatrix localPrec = pairPrecisions.get(i);
        
        RealVariable current = RealVariable.real(init[i+1]);
        variables.add(current);
        
        List<RealVariable> variables = new ArrayList<>();
        variables.add(prev);
        variables.add(current);
        
        CollisionFactor f = 
            useAnalytic ? 
            NormalFactor.newBinaryFactor(localPrec, current, prev) :
            NumericNormalFactor.withPrecision(variables, localPrec);
        
        result.add(f);
        prev = current;
      }
      
      return result;
    }
  }

}
