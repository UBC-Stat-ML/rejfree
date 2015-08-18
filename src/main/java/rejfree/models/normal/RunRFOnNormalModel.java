package rejfree.models.normal;

import org.jblas.DoubleMatrix;
import org.junit.Test;

import rejfree.local.LocalRFRunner;
import rejfree.models.normal.NormalChain.NormalChainModel;
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
  public LocalRFRunner runner = new LocalRFRunner();
  
  @Option
  public int nRepeats = 100;
  
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
      long seed = this.options.random.nextLong();
      DoubleMatrix exactSample = chain.exactSample();
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
      runner.init(modelSpec);
      runner.addMomentRayProcessor();
      runner.run();
      
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
