package rejfree;

import java.io.File;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TTest;
import org.jblas.DoubleMatrix;

import bayonet.rplot.PlotHistogram;

import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import rejfree.SimpleRFSampler.SimpleRFSamplerOptions;



public class TestBivariateNormDistribution
{

  public static void main(String[] args)
  {
    NormalEnergy energy = NormalEnergy.isotropic(2);
      
    SimpleRFSamplerOptions samplerOptions = new SimpleRFSamplerOptions();
    samplerOptions.collectRate = 0.0005;
    SimpleRFSampler sampler = SimpleRFSampler.initializeRFWithLBFGS(energy, samplerOptions);
    Random rand = new Random(134);
    int nIterations = 4000000;
    sampler.iterate(rand, nIterations);
    System.out.println("nSamples = " + sampler.getSamples().size());
    List<Double> 
      rfSamples = 
        //ref(nIterations, rand),
        normSqr(sampler.getSamples()),
      ref = ref(sampler.getSamples().size(), rand);
    
    PlotHistogram.from(rfSamples).toPDF(new File("hist-rejf.pdf"));
    PlotHistogram.from(ref).toPDF(new File("hist-truth.pdf"));
    
    MannWhitneyUTest test = new MannWhitneyUTest();
    System.out.println(test.mannWhitneyUTest(Doubles.toArray(rfSamples), Doubles.toArray(ref)));
    
    TTest test2 = new TTest();
    System.out.println(test2.tTest(Doubles.toArray(rfSamples), Doubles.toArray(ref)));
  }
  

  private static List<Double> ref(int nIterations, Random rand)
  {
    List<Double> result = Lists.newArrayList();
    for (int i = 0; i < nIterations; i++)
    {
      double normal1 = rand.nextGaussian();
      double normal2 = rand.nextGaussian();
      DoubleMatrix dm = new DoubleMatrix(new double[]{normal1,normal2});
      result.add(dm.dot(dm));
    }
    return result;
  }

  private static List<Double> normSqr(List<DoubleMatrix> samples)
  {
    List<Double> result = Lists.newArrayList();
    for (DoubleMatrix sample : samples)
      result.add(sample.dot(sample));
    return result;
  }

}
