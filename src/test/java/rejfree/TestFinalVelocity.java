package rejfree;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;

import bayonet.rplot.PlotHistogram;
import blang.validation.CheckStationarity;
import briefj.run.Mains;
import briefj.run.Results;
import rejfree.GlobalRFSampler.RFSamplerOptions;



public class TestFinalVelocity implements Runnable
{

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new TestFinalVelocity());
  }
  
  public void run()
  {
    for (double fixedTime : new double[]{0.0, 10})
    {
      final int dim = 2;
      final int nSamples = 100000;
      List<Double> betas = new ArrayList<Double>();
      List<Double> xs = new ArrayList<Double>();
      List<Double> norms = new ArrayList<Double>();
      DoubleMatrix
        outerProdFromExact = new DoubleMatrix(dim,dim),
        outerProdFromRF = new DoubleMatrix(dim,dim);
      
      RFSamplerOptions options = new RFSamplerOptions();
      Random rand = new Random(1);
      for (int i = 0; i < nSamples; i++)
      {
        NormalEnergy e = NormalEnergy.isotropic(dim);
        DoubleMatrix initPos = DoubleMatrix.randn(dim);
        
        outerProdFromExact.addi(initPos.mmul(initPos.transpose()));
        
        options.refreshRate = 0.0;
        options.collectRate = 0.0;
        GlobalRFSampler sampler = new GlobalRFSampler(e, initPos.dup() , options);
        
        sampler.iterate(rand , (int)fixedTime * 2+2);
  //      System.out.println(sampler.getTrajectory());
        List<DoubleMatrix> trajectory = sampler.getTrajectory();
        
  
        Pair<DoubleMatrix, DoubleMatrix> findFinal = findFinal(trajectory, fixedTime);
        DoubleMatrix x = findFinal.getLeft(),
          v = findFinal.getRight();
        double beta = Math.acos(x.dot(v) / x.norm2());
        betas.add(beta);
        xs.add(x.get(0));
        norms.add(x.norm2());
        
        outerProdFromRF.addi(x.mmul(x.transpose()));
      }
      outerProdFromExact.divi(nSamples);
      outerProdFromRF.divi(nSamples);
      
      System.out.println("Empirical covar from exact samples");
      System.out.println(outerProdFromExact);
      
      System.out.println("Empirical covar from RF samples");
      System.out.println(outerProdFromRF);
      
      PlotHistogram.from(norms).toPDF(Results.getFileInResultFolder("norm-trajLen=" + fixedTime + ".pdf"));
      PlotHistogram.from(betas).toPDF(Results.getFileInResultFolder("beta-trajLen=" + fixedTime + ".pdf"));
      PlotHistogram.from(xs).toPDF(Results.getFileInResultFolder("marg-trajLen=" + fixedTime + ".pdf"));
      
      System.out.println(nColls);
      
      List<Double> reference = new ArrayList<Double>();
      List<Double> refNorms = new ArrayList<Double>();
      for (int i = 0; i < nSamples; i++)
      {
        reference.add(rand.nextGaussian());
        DoubleMatrix mvn = DoubleMatrix.randn(dim);
        refNorms.add(mvn.norm2());
      }
      
      System.out.println("Marginal");
      System.out.println("p = " + CheckStationarity.tTest.pValue(xs, reference));
      System.out.println("p = " + CheckStationarity.higherMomentTTest.pValue(xs, reference));
      System.out.println("p = " + CheckStationarity.mannWhitneyTest.pValue(xs, reference));
      
      System.out.println("Norm");
      System.out.println("p = " + CheckStationarity.tTest.pValue(norms, refNorms));
      System.out.println("p = " + CheckStationarity.higherMomentTTest.pValue(norms, refNorms));
      System.out.println("p = " + CheckStationarity.mannWhitneyTest.pValue(xs, reference));
    }
  }
  
  static SummaryStatistics nColls = new SummaryStatistics();
  
  // x, v
  public static Pair<DoubleMatrix, DoubleMatrix> findFinal(List<DoubleMatrix> trajectory, double fixedTime)
  {
    double currentTimeConsumed = 0.0;
    for (int i = 0; i < trajectory.size() - 1; i++)
    {
      DoubleMatrix cur = trajectory.get(i),
           nxt = trajectory.get(i+1);
      
      DoubleMatrix diff = nxt.sub(cur);
      double norm = diff.norm2();
      
      if (currentTimeConsumed + norm > fixedTime)
      {
        nColls.addValue(i);
        
        double timeLeft = fixedTime - currentTimeConsumed;
        if (timeLeft < 0)
          throw new RuntimeException();
        
        DoubleMatrix velo = diff.div(diff.norm2());
        DoubleMatrix point = cur.add(velo.mul(timeLeft));
        return Pair.of(point, velo);
      }
      
      currentTimeConsumed += norm;
    }
    throw new RuntimeException();
  }

}
