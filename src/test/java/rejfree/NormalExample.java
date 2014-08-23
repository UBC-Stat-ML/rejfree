package rejfree;

import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;

import com.google.common.collect.Lists;

import bayonet.rplot.PlotHistogram;
import briefj.run.Mains;
import briefj.run.Results;





public class NormalExample implements Runnable
{
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new NormalExample());
  }

  @Override
  public void run()
  {
    DoubleMatrix covar = new DoubleMatrix(new double[][]{{1,0},{0,10}});
    NormalEnergy energy = NormalEnergy.
      withCovariance(covar); 
      //isotropic(2);
    SimpleRFSampler sampler = new SimpleRFSampler(energy);
    Random rand = new Random(134);
    sampler.iterate(rand, 1000);
    PlotTrajectory pt = new PlotTrajectory(sampler.getTrajectory(), 0, 1);
    pt.setWorkFolder(Results.getFolderInResultFolder("r-scripts"));
    pt.toPDF(Results.getFileInResultFolder("trajectory.pdf"));
    System.out.println(sampler.getCollisionToRefreshmentRatio().getMean());
    System.out.println(sampler.getCollectedPerEvent().getMean());
    
    for (int d = 0; d < 2; d++)
    {
      List<Double> coordinates = Lists.newArrayList();
      for (DoubleMatrix sample : sampler.getSamples())
        coordinates.add(sample.get(d));
      PlotHistogram.from(coordinates).toPDF(Results.getFileInResultFolder("hist-" + d + ".pdf"));
    }
  }
}
