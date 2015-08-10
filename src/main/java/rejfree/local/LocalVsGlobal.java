package rejfree.local;

import java.util.concurrent.TimeUnit;

import org.jblas.DoubleMatrix;
import org.junit.Test;

import com.google.common.base.Stopwatch;

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
  public long maxRunningTimeMilli = Long.MAX_VALUE;
  
  @Option
  public int maxSteps = 1000;
  
  private NormalChain chain;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new LocalVsGlobal());
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
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data, isLocal);
      ProbabilityModel model = new ProbabilityModel(modelSpec);
      LocalRFSampler local = new LocalRFSampler(model, rfOptions);
      Stopwatch watch = Stopwatch.createStarted();
      local.iterate(this.options.random, maxSteps, Double.POSITIVE_INFINITY, maxRunningTimeMilli);
      long elapsed = watch.elapsed(TimeUnit.MILLISECONDS);
      output.printWrite("time", "seed", seed, "timeMilli", elapsed, "nCollisions", local.getNCollisions());
      
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable variable = modelSpec.variables.get(d);
        double truth = chain.covarMatrix.get(d, d);
        double estimate = local.getVarEstimate(variable);
        double error = Math.abs(truth - estimate);
        output.printWrite("results", "dim", d, "seed", seed, "absError", error, "relError", (error/truth), "truth", truth, "estimate", estimate);
      }
    }
    
    output.close();
  }

}
