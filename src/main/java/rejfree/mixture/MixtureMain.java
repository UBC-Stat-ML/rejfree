package rejfree.mixture;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;

import rejfree.PlotTrajectory;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.local.LocalRFSampler;
import rejfree.local.TrajectoryRay;
import rejfree.mixture.MixtureModel.MixtureSpec;
import rejfree.processors.RayProcessor;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;


/**
 * Warning: this was not implemented with the goal of efficiency or
 * generality in mind, but rather to quickly create intuitive plots
 * explaining the execution of the algorithm.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class MixtureMain implements Runnable
{
  @OptionSet(name = "mix")
  public final MixtureModel model = new MixtureModel();
  
  @OptionSet(name = "rf")
  public final LocalRFRunnerOptions localRFRunnerOption = new LocalRFRunnerOptions();

  @Override
  public void run()
  {
    OutputManager output = Results.getGlobalOutputManager();
    
    model.densityPlot().toPDF(Results.getFileInResultFolder("targetDensity.pdf"));
    
    MixtureSpec modelSpec = model.getModelSpec();
    LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);
    rfRunner.init(modelSpec);
    List<DoubleMatrix> collisionPositions = new ArrayList<>();
    rfRunner.sampler.addRayProcessor(new RayProcessor() {
      
      @Override
      public void processRay(RealVariable var, TrajectoryRay ray, double time,
          LocalRFSampler sampler)
      {
        sampler.updateAllVariables(time);
        DoubleMatrix current = new DoubleMatrix(2);
        current.put(0, modelSpec.meanForComponentMean0.getValue());
        current.put(1, modelSpec.meanForComponentMean1.getValue());
        collisionPositions.add(current);
      }
      
      @Override
      public void init(LocalRFSampler sampler)
      {
      }
    });
    rfRunner.run();
    
    PlotTrajectory plot = PlotTrajectory.fromFirstTwoDimensions(collisionPositions);
    plot.xBounds = new double[]{-2,2};
    plot.yBounds = new double[]{-2,2};
    plot.toPDF(Results.getFileInResultFolder("traj.pdf"));
    
    output.close();
  }
  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new MixtureMain());
  }
  
}
