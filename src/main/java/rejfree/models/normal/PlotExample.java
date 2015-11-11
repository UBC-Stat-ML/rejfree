package rejfree.models.normal;

import java.util.List;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.jblas.DoubleMatrix;

import rejfree.PlotTrajectory;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.mixture.MixtureMain.LocAtCollisionsProcessor;
import bayonet.rplot.PlotContour;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class PlotExample implements Runnable
{
  @OptionSet(name = "rf")
  public LocalRFRunnerOptions localRFRunnerOption = new LocalRFRunnerOptions();
  
  @Option
  public int bound = 3;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new PlotExample());
  }
  
  public class Spec
  {
    RealVariable variable0 = new RealVariable(1), variable1 = new RealVariable(1);
    
    @DefineFactor
    public NormalFactor factor = NormalFactor.newBinaryFactor(precision(), variable0, variable1);

    private DoubleMatrix precision()
    {
      DoubleMatrix result = new DoubleMatrix(2,2);
      result.put(0,0,1);
      result.put(1,1,1);
      return result;
    }
  }
  
  public PlotContour densityPlot()
  {
    Spec modelSpec = new Spec();
    ProbabilityModel m = new ProbabilityModel(modelSpec);
    
    PlotContour plot = new PlotContour(new MultivariateFunction() {
      
      @Override
      public double value(double[] point)
      {
        modelSpec.variable0.setValue(point[0]);
        modelSpec.variable1.setValue(point[1]);
        return m.logDensity();
      }
    });
    final double 
      min = -bound,
      max = bound;
    
    plot.min_x = min;
    plot.min_y = min;
    plot.max_x = max;
    plot.max_y = max;
    
    return plot;
  }


  @Override
  public void run()
  {
    PlotContour densityPlot = densityPlot();
    densityPlot.toPDF(Results.getFileInResultFolder("density.pdf"));
    
    Spec modelSpec = new Spec();
    LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);
    rfRunner.init(modelSpec);
    
    LocAtCollisionsProcessor lcp = new LocAtCollisionsProcessor(modelSpec.variable0, modelSpec.variable1);
    rfRunner.sampler.addRayProcessor(lcp);
    rfRunner.run();
  
//    int start = findGoodIndex(collisionPositions);
    List<DoubleMatrix> toPlot =
//        collisionPositions.subList(start-2, start+4)
        lcp.locationsAtCollisions;
    
    PlotTrajectory plot = PlotTrajectory.fromFirstTwoDimensions(toPlot);
    plot.xBounds = new double[]{-bound,bound};
    plot.yBounds = new double[]{-bound,bound};
    plot.toPDF(Results.getFileInResultFolder("traj.pdf"));
    
  }
}
