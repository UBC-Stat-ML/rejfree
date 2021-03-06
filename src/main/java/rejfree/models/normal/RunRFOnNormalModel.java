package rejfree.models.normal;


import java.util.List;

import org.jblas.DoubleMatrix;
import org.junit.Test;

import rejfree.PlotTrajectory;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.models.normal.NormalChain.NormalChainModel;
import rejfree.processors.LocAtCollisionsProcessor;
import bayonet.coda.EffectiveSize;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class RunRFOnNormalModel implements Runnable
{
  @OptionSet(name = "modelOptions")
  public NormalChainOptions options = new NormalChainOptions();
  
  @OptionSet(name = "rfRunner")
  public LocalRFRunnerOptions runnerOptions = new LocalRFRunnerOptions();
  
  @Option
  public int nRepeats = 100;
  
  @Option
  public int nColToPlot = 0;
  
  @Option
  public double plotBound = 3;
  
  private NormalChain chain;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new RunRFOnNormalModel());
  }
  
  @Test
  public void run()
  {
    org.jblas.util.Random.seed(options.random.nextLong());
    chain = new NormalChain(options);
    OutputManager output = new OutputManager();
    output.setOutputFolder(Results.getResultFolder());
    
    for (int rep = 0; rep < nRepeats; rep++)
    {
      LocalRFRunner runner = new LocalRFRunner(runnerOptions);
      long seed = this.options.random.nextLong();
      DoubleMatrix exactSample = chain.exactSample();
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
      runner.init(modelSpec);
      runner.addMomentRayProcessor();
      runner.addSaveAllRaysProcessor();
      
      LocAtCollisionsProcessor lcp = new LocAtCollisionsProcessor(modelSpec.variables.get(0), modelSpec.variables.get(1));
      if (nColToPlot > 0)
        runner.sampler.addRayProcessor(lcp);
      
      runner.run();
      
      if (nColToPlot > 0)
      {
        List<DoubleMatrix> list = lcp.locationsAtCollisions;
        PlotTrajectory plot = PlotTrajectory.fromFirstTwoDimensions(list.subList(0, Math.min(list.size(), nColToPlot)));
        plot.xBounds = new double[]{-plotBound,plotBound};
        plot.yBounds = new double[]{-plotBound,plotBound};
        plot.toPDF(Results.getFileInResultFolder("traj-" + rep + ".pdf"));
      }
      
      RealVariable aVar = modelSpec.variables.get(3);
      double rfESS = EffectiveSize.effectiveSize(runner.saveRaysProcessor.convertToSample(aVar, 4.0));
      Results.getGlobalOutputManager().printWrite("essPerSec", "value", rfESS);
      
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable variable = modelSpec.variables.get(d);
        double truth = chain.covarMatrix.get(d, d);
        double estimate = runner.momentRayProcessor.getVarianceEstimate(variable);
        double error = Math.abs(truth - estimate);
        output.printWrite("results", 
            "dim", d, 
            "rep", rep,
            "seed", seed, 
            "absError", error, 
            "relError", (error/truth), 
            "truth", truth, 
            "estimate", estimate);
      }
    }
    
    output.close();
  }

}
