package rejfree.local;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFSampler;
import bayonet.math.JBlasUtils;
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
  
  DoubleMatrix outerProduct = new DoubleMatrix(3,3);
  
  public static void main(String [] args)
  {
    new TestLocalOnNormalDist().test();
  }
  
  @Test
  public void test()
  {
    for (boolean useLocalRef : new boolean[]{true,false})
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
      
      
      RFSamplerOptions options = new RFSamplerOptions();
      options.useLocalRefreshment = useLocalRef;
      options.refreshRate = 0.0001;
      options.collectRate = 10.0;
      LocalRFSampler local = new LocalRFSampler(model, options);
      
      local.processors.clear();
      local.processors.add(this);
      Random rand = new Random(1);
      local.iterate(rand , 200000, Double.POSITIVE_INFINITY);
      
      for (int i = 0; i < 3; i++)
      {
        System.out.println(stats.get(i).getVariance());
        System.out.println(stats.get(i).getMean());
        Assert.assertEquals(fullCovar.get(i,i), stats.get(i).getVariance(), 0.02);
        Assert.assertEquals(0.0, stats.get(i).getMean(), 0.02);
        System.out.println("---");
      }
      
      outerProduct.divi(stats.get(0).getN());
      System.out.println("Empirical covar matrix");
      System.out.println(outerProduct);
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
    NumericNormalFactor f12 = NumericNormalFactor.withCovariance(Lists.newArrayList(v1, v2), localCovar);
    
    @DefineFactor
    NumericNormalFactor f23 = NumericNormalFactor.withCovariance(Lists.newArrayList(v2, v3), localCovar);
  }

  private List<SummaryStatistics> stats = new ArrayList<>();
  @Override
  public void process(ProcessorContext context)
  {
    DoubleMatrix cur = new DoubleMatrix(3);
    for (int i = 0; i < 3; i++)
    {
      final double value = modelSpec.vs[i].getValue();
      stats.get(i).addValue(value);
      cur.put(i, value);
    }
    outerProduct.addi(cur.mmul(cur.transpose()));
    
  }
}
