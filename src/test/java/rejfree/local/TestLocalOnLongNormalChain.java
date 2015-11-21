package rejfree.local;

import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import rejfree.RFSamplerOptions;
import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.local.LocalRFSampler;
import rejfree.models.normal.NormalChain;
import rejfree.models.normal.NormalChainOptions;
import rejfree.models.normal.NormalChain.NormalChainModel;
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
  
  public void testMCAvgs(RefreshmentMethod refreshMethod)
  {
    System.out.println("Checking MC averages");
    LocalRFRunner runner = new LocalRFRunner();
    runner.options.rfOptions.refreshRate = 2.0;
    runner.options.rfOptions.collectRate = 1.0;
    runner.options.rfOptions.refreshmentMethod = refreshMethod;
    
    final DoubleMatrix 
      outerSumsRF = new DoubleMatrix(chain.dim(),chain.dim());
    int [] i = new int[1];
    DoubleMatrix exactSample = chain.exactSample();
    NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
    runner.init(modelSpec);
    runner.addMomentRayProcessor();
    runner.addSaveAllRaysProcessor();
    runner.sampler.addPointProcessor(new Processor()
    {
      @Override
      public void process(ProcessorContext context)
      {
        DoubleMatrix x = currentPositionToVector(modelSpec);
        outerSumsRF.addi(x.mmul(x.transpose()));
        i[0]++;
      }
    });
    runner.sampler.iterate(this.options.random, 1_000_000, Double.POSITIVE_INFINITY);
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
      double meanEstimate = runner.momentRayProcessor.getMeanEstimate(variable);
      double varEstimate = runner.momentRayProcessor.getVarianceEstimate(variable);
      System.out.println(meanEstimate + "\t" + varEstimate);
      Assert.assertTrue(meanEstimate < 0.05);
      Assert.assertTrue(Math.abs(varEstimate - chain.covarMatrix.get(d,d)) < 0.05);
    } 
    
    System.out.println("from regular spaced samples:");
    final double delta = 1.0;
    for (int d = 0; d < chain.dim(); d++)
    {
      RealVariable variable = modelSpec.variables.get(d);
      List<Double> converted = runner.saveRaysProcessor.convertToSample(variable, delta);
      SummaryStatistics stats = new SummaryStatistics();
      for (double cur : converted)
        stats.addValue(cur);
      
      double meanEstimate = stats.getMean();
      double varEstimate = stats.getVariance();
      System.out.println(meanEstimate + "\t" + varEstimate);
      Assert.assertTrue(meanEstimate < 0.05);
      Assert.assertTrue(Math.abs(varEstimate - chain.covarMatrix.get(d,d)) < 0.05);
    } 
  }
  
  public void testInvariance(RefreshmentMethod method)
  {
    System.out.println("Checking invariance");
    RFSamplerOptions options = new RFSamplerOptions();
    options.refreshRate = 1.0;
    options.collectRate = 0.0;
    options.refreshmentMethod = method;

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
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
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
    options.nPairs = 3;
    
    chain = new NormalChain(options);
    
    System.out.println("true covariance matrix:\n" + chain.covarMatrix);

    for (RefreshmentMethod method : RefreshmentMethod.values())
    {
      System.out.println("=====");
      System.out.println("refreshmentMethod = " + method);
      testMCAvgs(method);
      testInvariance(method);
    }
  }

  public static void main(String[] args)
  {
    new TestLocalOnLongNormalChain().run();
  }


}
