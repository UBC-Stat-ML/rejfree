package rejfree.local;

import org.jblas.DoubleMatrix;
import org.junit.Test;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.NormalChain.NormalChainModel;
import blang.ProbabilityModel;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class LocalVsGlobal implements Runnable
{
  @OptionSet(name = "modelOptions")
  public NormalChainOptions options = new NormalChainOptions();
  
  @OptionSet(name = "rfOptions")
  public RFSamplerOptions rfOptions = new RFSamplerOptions();
  
  @Option
  public int nRepeats = 100;
  
  @Option
  public boolean isLocal = false;
  
  @Option
  public double timeInterval = 1000.0;
  
  private NormalChain chain;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new LocalVsGlobal());
  }

  @Test
  public void run()
  {
    org.jblas.util.Random.seed(options.random.nextLong()); // not needed?
    chain = new NormalChain(options);
    OutputManager output = new OutputManager();
    output.setOutputFolder(Results.getResultFolder());
    
    for (int rep = 0; rep < nRepeats; rep++)
    {
      long seed = this.options.random.nextLong();
      DoubleMatrix exactSample = chain.exactSample();
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data, isLocal);
      ProbabilityModel model = new ProbabilityModel(modelSpec);
      LocalRFSampler local = new LocalRFSampler(model, rfOptions);
      local.iterate(this.options.random, Integer.MAX_VALUE, timeInterval);
      
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable var = modelSpec.variables.get(d);
        double trueVar = chain.covarMatrix.get(d, d);
        double estimate = local.getVarEstimate(var);
        double error = Math.abs(trueVar - estimate);
        output.printWrite("results", "dim", d, "seed", seed, "varError", error);
      }
    }
    
    output.close();
  }

}
