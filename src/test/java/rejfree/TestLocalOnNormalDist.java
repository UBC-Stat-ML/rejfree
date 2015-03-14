package rejfree;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;

import rejfree.SimpleRFSampler.SimpleRFSamplerOptions;
import bayonet.math.JBlasUtils;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.RealVariable;



public class TestLocalOnNormalDist implements Processor
{
  private static DoubleMatrix localCovar = new DoubleMatrix(new double[][]{{1,0.5},{0.5,1}});
  private ProbabilityModel model;
  private ChainModel modelSpec;
  
  public static void main(String [] args)
  {
    new TestLocalOnNormalDist().test();
  }
  
  @Test
  public void test()
  {
    modelSpec = new ChainModel();
    model = new ProbabilityModel(modelSpec);
    
    System.out.println(model);
    
    DoubleMatrix localPrecision = JBlasUtils.inversePositiveMatrix(localCovar);
    System.out.println(localPrecision);
    
    DoubleMatrix fullPrecision = new DoubleMatrix(3, 3);
    for (int i = 0; i < 2; i++)
    {
      for (int row = 0; row < localPrecision.rows; row++)
        for (int col = 0; col < localPrecision.columns; col++)
          fullPrecision.put(row + i, col + i, fullPrecision.get(row + i, col + i) + localPrecision.get(row, col));
    }
    
    System.out.println(fullPrecision);
    DoubleMatrix fullCovar = JBlasUtils.inversePositiveMatrix(fullPrecision);
    System.out.println(fullCovar);
    
    for (int i = 0; i < 3; i++)
      stats.add(new SummaryStatistics());
    
    
//    // standard
//    MCMCFactory factory = new MCMCFactory();
//    factory.mcmcOptions.nMCMCSweeps = 100000;
//    MCMCAlgorithm mcmc = factory.build(model);
//    mcmc.run();
    
    SimpleRFSamplerOptions options = new SimpleRFSamplerOptions();
    options.refreshRate = 0.0;
    LocalRFSampler local = new LocalRFSampler(model, options);
    
    
    local.processors.add(this);
    Random rand = new Random(1);
    local.iterate(rand , 100000);
    
    for (int i = 0; i < 3; i++)
    {
      System.out.println(stats.get(i).getVariance());
      System.out.println(stats.get(i).getMean());
      Assert.assertEquals(stats.get(i).getVariance(), fullCovar.get(i,i), 0.01);
      Assert.assertEquals(stats.get(i).getMean(), 0.0, 0.01);
      System.out.println("---");
    }
    

  }
  
  
  public static class ChainModel
  {
    RealVariable 
      v1 = RealVariable.real(),
      v2 = RealVariable.real(),
      v3 = RealVariable.real();
    
    public RealVariable [] vs = new RealVariable[]{v1, v2, v3};

    @DefineFactor
    NormalFactor f12 = NormalFactor.withCovariance(Lists.newArrayList(v1, v2), localCovar);
    
    @DefineFactor
    NormalFactor f23 = NormalFactor.withCovariance(Lists.newArrayList(v2, v3), localCovar);
  }


  private List<SummaryStatistics> stats = new ArrayList<>();
  @Override
  public void process(ProcessorContext context)
  {
    for (int i = 0; i < 3; i++)
      stats.get(i).addValue(modelSpec.vs[i].getValue());
  }
}
