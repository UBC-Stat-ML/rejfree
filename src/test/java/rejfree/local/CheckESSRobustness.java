package rejfree.local;

import java.util.Collections;

import org.jblas.DoubleMatrix;
//import org.junit.Test;

import rejfree.models.normal.NormalChain;
import rejfree.models.normal.NormalChainOptions;
import rejfree.models.normal.NormalChain.NormalChainModel;
import bayonet.coda.EffectiveSize;
import blang.variables.RealVariable;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;




public class CheckESSRobustness implements Runnable
{
  @OptionSet(name = "modelOptions")
  public NormalChainOptions options = new NormalChainOptions();
  
  @OptionSet(name = "rfRunner")
  public LocalRFRunner runner = new LocalRFRunner();
  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new CheckESSRobustness());
  }
  
//  @Test
  public void run()
  {
    org.jblas.util.Random.seed(options.random.nextLong());
    NormalChain chain = new NormalChain(options);
    
    DoubleMatrix exactSample = chain.exactSample();
    NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
    runner.init(modelSpec);
    
    RealVariable aVar = modelSpec.variables.get(0);
    runner.addSaveRaysProcessor(Collections.singleton(aVar));
    runner.addSaveSamplesProcessor(Collections.singletonList(aVar));
    runner.run();
    
    double sampledESS = EffectiveSize.effectiveSize(runner.saveSamplesProcessor.samples.get(0));
    Results.getGlobalOutputManager().printWrite("ess", "sampled", true, "delta", (1.0/runner.options.rfOptions.collectRate), "ess", sampledESS);
    
    for (double delta = 1024; delta > 0.25; delta /= 2)
    {
      double ess = EffectiveSize.effectiveSize(runner.saveRaysProcessor.convertToSample(aVar, delta));
      Results.getGlobalOutputManager().printWrite("ess", "sampled", false, "delta", delta, "ess", ess);
    }
    
    Results.getGlobalOutputManager().flush();
    
  }

}
