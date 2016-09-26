package rejfree;

import java.util.ArrayList;
import java.util.List;

import bayonet.rplot.PlotHistogram;
import blang.variables.RealVariable;
import briefj.run.Mains;
import briefj.run.Results;
import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.local.LocalRFSampler;
import rejfree.local.TrajectoryRay;
import rejfree.processors.RayProcessor;
import rejfree.scalings.EstimateESSByDimensionality;

public class TrajLenScaling implements Runnable
{

  public static void main(String [] args) 
  {
    Mains.instrumentedRun(args, new TrajLenScaling());
  }
  
  public void run() 
  {
    for (boolean sparse : new boolean[]{true,false}) 
    {
      LocalRFRunnerOptions rOptions = new LocalRFRunnerOptions();
      rOptions.maxSteps = Integer.MAX_VALUE;
      rOptions.maxTrajectoryLength = 100.0;
      rOptions.rfOptions.collectRate = 0.0;
      rOptions.rfOptions.refreshmentMethod = RefreshmentMethod.GLOBAL;
      for (int index = 0; index < 15; index++)
      {
        int dim = 1 << index;
        LocalRFRunner runner = new LocalRFRunner(rOptions);
        runner.init(new EstimateESSByDimensionality.ModelSpec(dim, sparse));
        SegmentLengthProcessor sp = new SegmentLengthProcessor();
        runner.sampler.rayProcessors.add(sp);
        runner.run();
        PlotHistogram.from(sp.deltas).toPDF(Results.getFileInResultFolder("sparse=" + sparse + ",dim=" + dim + ".pdf"));
      }
    }
  }
  
  public static class SegmentLengthProcessor implements RayProcessor
  {
    double lastEvent = 0.0;
    List<Double> deltas = new ArrayList<>();

    @Override
    public void init(LocalRFSampler sampler)
    {
    }

    @Override
    public void processRay(RealVariable var, TrajectoryRay ray, double time, LocalRFSampler sampler)
    {
      if (time == lastEvent)
        return;
      deltas.add(time - lastEvent);
      lastEvent = time;
    }
  }
}
