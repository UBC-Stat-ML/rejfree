package rejfree.models.normal;

import java.io.File;
import java.util.List;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.jblas.DoubleMatrix;

import rejfree.PlotTrajectory;
import rejfree.StaticUtils;
import rejfree.local.CollisionFactor;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.processors.LocAtCollisionsProcessor;
import bayonet.opt.DifferentiableFunction;
import bayonet.rplot.PlotContour;
import bayonet.rplot.RUtils;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.OutputManager;
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
    RealVariable variable0 = new RealVariable(1), variable1 = new RealVariable(0);
    
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
  
  public static boolean hackInitMode = false;

  @Override
  public void run()
  {
    hackInitMode = true;
    PlotContour densityPlot = densityPlot();
    densityPlot.toPDF(Results.getFileInResultFolder("density.pdf"));
    
    Spec modelSpec = new Spec();
    LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);
    rfRunner.init(modelSpec);
    
    LocAtCollisionsProcessor lcp = new LocAtCollisionsProcessor(modelSpec.variable0, modelSpec.variable1);
    rfRunner.sampler.addRayProcessor(lcp);
    rfRunner.run();
  
    List<DoubleMatrix> toPlot =
        lcp.locationsAtCollisions;
    
    PlotTrajectory plot = PlotTrajectory.fromFirstTwoDimensions(toPlot);
    plot.xBounds = new double[]{-bound,bound};
    plot.yBounds = new double[]{-bound,bound};
    plot.toPDF(Results.getFileInResultFolder("traj.pdf"));
    
    for (int i = 0; i < lcp.locationsAtCollisions.size() - 1; i++)
    {
      DoubleMatrix 
        start = lcp.locationsAtCollisions.get(i),
        end = lcp.locationsAtCollisions.get(i+1);
      File dir = Results.getFolderInResultFolder("seg" + i);
      plotSegmentCurves(start, end, energyFunction(modelSpec.factor), dir, 1000);
    }
  }
  
  public static DifferentiableFunction energyFunction(CollisionFactor factor)
  {
    return new DifferentiableFunction() {
      
      private void setTo(double [] x)
      {
        for (int i = 0; i < factor.nVariables(); i++)
          factor.getVariable(i).setValue(x[i]);
      }
      
      @Override
      public double valueAt(double[] x)
      {
        setTo(x);
        return -factor.logDensity();
      }
      
      @Override
      public int dimension()
      {
        return factor.nVariables();
      }
      
      @Override
      public double[] derivativeAt(double[] x)
      {
        setTo(x);
        DoubleMatrix gradient = factor.gradient();
        gradient.muli(-1);
        return gradient.data;
      }
    };
  }
  
  public static void plotSegmentCurves(DoubleMatrix start, DoubleMatrix end, DifferentiableFunction energyFunction, File dir, int nPoints)
  {
    OutputManager out = new OutputManager();
    out.setOutputFolder(dir);
    
    // plot the energy
    DoubleMatrix diff = end.sub(start);
    double time = diff.norm2();
    DoubleMatrix direction = diff.div(time);
    out.write("collisionTime", "value", time);
    
    double deltaT = 1.5*time/nPoints;
    for (int i = 0; i < nPoints; i++)
    {
      double curT = deltaT*i;
      DoubleMatrix position = StaticUtils.position(start, direction, curT);
      
      double energy = energyFunction.valueAt(position.data);
      
      double dotProd = direction.dot(new DoubleMatrix(energyFunction.derivativeAt(position.data)));
      
      out.write("energy", "t", curT, "name", "energy",     "value", energy);
      out.write("deriv", "t", curT, "name", "derivative", "value", dotProd);
      out.write("deriv", "t", curT, "name", "intensity",  "value", Math.max(0,dotProd));
    }
    out.close();
    
    RUtils.callGeneratedRScript("/rejfree/plotEnergyAndIntensity.txt", dir);
  }
}
