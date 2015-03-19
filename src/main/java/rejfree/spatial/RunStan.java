package rejfree.spatial;


import java.io.File;

import bayonet.bugs.StanWrapper;
import binc.Command;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class RunStan implements Runnable
{
  @Option
  public File stanHome = new File("/Users/bouchard/Documents/workspace-cpp/cmdstan");
  
  @Option
  public File data = new File("/Users/bouchard/Documents/experiments/rejfree-spatial/results/all/2015-03-17-20-33-01-Q6rfxVgn.exec/output.r");

  @Option
  public int thin = 1;

  @Option
  public int nSamples = 1000;

  @Option
  public int nWarmUp = 1000;

  @Option
  public int seed = 1;
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new RunStan());
  }

  @Override
  public void run()
  {
    String stanModel = BriefIO.resourceToString("/rejfree/spatialModel.stan");
    File compiled = StanWrapper.compile(stanModel, stanHome);
    Command.byPath(compiled)
      .ranIn(Results.getResultFolder())
      .withStandardOutMirroring()
      .withArgs("sample "
          + "thin=" + thin  + " "
          + "num_samples=" + nSamples + " "
          + "num_warmup=" + nWarmUp + " "
          + "random seed=" + seed  + " "
          + "data file=" + data.getAbsolutePath())
      .call();
  }
  
}
