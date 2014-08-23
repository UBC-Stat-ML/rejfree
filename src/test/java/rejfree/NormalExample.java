package rejfree;

import java.util.Random;

import org.jblas.DoubleMatrix;

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
    sampler.setCurrentPosition(new DoubleMatrix(new double[]{1,2}));
    Random rand = new Random(134);
    sampler.iterate(rand, 1000);
    System.out.println(sampler.getTrajectory());
    PlotTrajectory pt = new PlotTrajectory(sampler.getTrajectory(), 0, 1);
    pt.setWorkFolder(Results.getFolderInResultFolder("r-scripts"));
    pt.toPDF(Results.getFileInResultFolder("trajectory.pdf"));
  }
}
