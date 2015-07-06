package rejfree.local;

import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFSampler;
import rejfree.local.NormalChain.NormalChainModel;
import blang.ProbabilityModel;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.RealVariable;
import briefj.opt.Option;



public class TestLocalOnLongNormalChain implements Runnable
{
  @Option
  public NormalChainOptions options = new NormalChainOptions();
  
  private NormalChain chain;

  public static double maxError(DoubleMatrix guess, DoubleMatrix truth)
  {
    return (guess.sub(truth)).normmax();
  }
  
  public void testMCAvgs()
  {
    RFSamplerOptions options = new RFSamplerOptions();
    options.refreshRate = 2.0;
    options.collectRate = 1.0;
    final DoubleMatrix 
      outerSumsRF = new DoubleMatrix(chain.dim(),chain.dim());
    int [] i = new int[1];
    DoubleMatrix exactSample = chain.exactSample();
    NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data, false);
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
    local.iterate(this.options.random, 1000000, Double.POSITIVE_INFINITY);
    outerSumsRF.divi(i[0]);
    System.out.println("empirical from MC avg (nSamples = " + i[0] + ")");
    System.out.println(outerSumsRF);
    
    double maxErr = maxError(chain.covarMatrix, outerSumsRF);
    System.out.println("maxError = " + maxErr);
    Assert.assertTrue(maxErr < 0.05);
    
    System.out.println("from line integral estimator:");
    for (int d = 0; d < chain.dim(); d++)
    {
      RealVariable variable = modelSpec.variables.get(d);
      double meanEstimate = local.getMeanEstimate(variable);
      double varEstimate = local.getVarEstimate(variable);
      System.out.println(meanEstimate + "\t" + varEstimate);
      Assert.assertTrue(meanEstimate < 0.05);
      Assert.assertTrue(Math.abs(varEstimate - chain.covarMatrix.get(d,d)) < 0.05);
    } 
  }
  
  public void testInvariance()
  {
    RFSamplerOptions options = new RFSamplerOptions();
    options.refreshRate = 1.0;
    options.collectRate = 0.0;

    double fixedTime = 10;
    int nRepeats = 20000;
    
    DoubleMatrix 
      outerSumsExact = new DoubleMatrix(chain.dim(),chain.dim()),
      outerSumsRF = new DoubleMatrix(chain.dim(),chain.dim());
    
    for (int i = 0; i < nRepeats; i++)
    {
      DoubleMatrix exactSample = chain.exactSample();
      outerSumsExact.addi(exactSample.mmul(exactSample.transpose()));
      
      DoubleMatrix rfSampl = null;
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data, true);
      ProbabilityModel model = new ProbabilityModel(modelSpec);
      LocalRFSampler local = new LocalRFSampler(model, options);
      local.iterate(this.options.random, Integer.MAX_VALUE, fixedTime);
      rfSampl = currentPositionToVector(modelSpec);
      outerSumsRF.addi(rfSampl.mmul(rfSampl.transpose()));
    }
    outerSumsExact.divi(nRepeats);
    outerSumsRF.divi(nRepeats);
    
    System.out.println("Empirical covar from exact samples");
    System.out.println(outerSumsExact);
    
    System.out.println("Empirical covar from RF samples");
    System.out.println(outerSumsRF);
    
    double maxErr = maxError(chain.covarMatrix, outerSumsRF);
    System.out.println("maxError = " + maxErr);
    Assert.assertTrue(maxErr < 0.05);
  }
  
  static DoubleMatrix currentPositionToVector(NormalChainModel model)
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
    org.jblas.util.Random.seed(options.random.nextLong());
    
    chain = new NormalChain(options);
    testMCAvgs();
    testInvariance();
  }

  public static void main(String[] args)
  {
    new TestLocalOnLongNormalChain().run();
  }


}
