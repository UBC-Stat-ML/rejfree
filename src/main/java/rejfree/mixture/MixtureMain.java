package rejfree.mixture;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashSet;
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
  
  public static class LocAtCollisionsProcessor implements RayProcessor
  {
    public final LinkedHashSet<RealVariable> variables;
    public final List<DoubleMatrix> locationsAtCollisions = new ArrayList<>();
    private double prevTime = -1;

    public LocAtCollisionsProcessor(RealVariable ... variables)
    {
      this.variables = new LinkedHashSet<>(Arrays.asList(variables));
    }
    
    public LocAtCollisionsProcessor(Collection<RealVariable> variables)
    {
      this.variables = new LinkedHashSet<RealVariable>(variables);
    }

    @Override
    public void init(LocalRFSampler sampler)
    {
    }

    @Override
    public void processRay(RealVariable var, TrajectoryRay ray, double time,
        LocalRFSampler sampler)
    {
      if (time == prevTime)
        return;
      
      if (!variables.contains(var))
        return;
      
      DoubleMatrix current = new DoubleMatrix(variables.size());
      int i = 0;
      for (RealVariable aVar : variables)
      {
        sampler.updateVariable(aVar, time);
        current.put(i, aVar.getValue());
        i++;
      }
      locationsAtCollisions.add(current);
      
      prevTime = time;
    }
    
  }

  @Override
  public void run()
  {
    OutputManager output = Results.getGlobalOutputManager();
    
    model.densityPlot().toPDF(Results.getFileInResultFolder("targetDensity.pdf"));
    
    MixtureSpec modelSpec = model.getModelSpec();
    LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);
    rfRunner.init(modelSpec);
    
    LocAtCollisionsProcessor lcp = new LocAtCollisionsProcessor(modelSpec.meanForComponentMean0, modelSpec.meanForComponentMean1);
    rfRunner.sampler.addRayProcessor(lcp);
    rfRunner.run();
  
//    int start = findGoodIndex(collisionPositions);
    List<DoubleMatrix> toPlot =
//        collisionPositions.subList(start-2, start+4)
        lcp.locationsAtCollisions;
    
    PlotTrajectory plot = PlotTrajectory.fromFirstTwoDimensions(toPlot);
    plot.xBounds = new double[]{-2,2};
    plot.yBounds = new double[]{-2,2};
    plot.toPDF(Results.getFileInResultFolder("traj.pdf"));
    
    output.close();
  }
  
  int findGoodIndex(List<DoubleMatrix> collisionPositions)
  {
    double max = Double.NEGATIVE_INFINITY;
    int argmax = -1;
    
    
    for (int i = 0; i < collisionPositions.size() - 1; i++)
    {
      DoubleMatrix
        cur = collisionPositions.get(i),
        nxt = collisionPositions.get(i+1);
      
      double distance = cur.distance2(nxt);
      
      if (distance > max)
      {
        max = distance;
        argmax = i;
      }
    }
    
    return argmax;
  }

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new MixtureMain());
  }
  
}
