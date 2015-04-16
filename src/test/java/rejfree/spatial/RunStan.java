package rejfree.spatial;


import java.io.File;

import bayonet.bugs.StanWrapper;
import binc.Command;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.run.Results;



public class RunStan 
{
  @Option
  public File stanHome;
  
  public final SpatialMainOptions mainOptions;

  @Option
  public int thin = 1;

  @Option
  public int nWarmUp = 1000;
  
  public RunStan(SpatialMainOptions mainOptions)
  {
    super();
    this.mainOptions = mainOptions;
  }

  public void run()
  {
    String stanModel = BriefIO.resourceToString("/rejfree/spatialModel.stan");
    File compiled = StanWrapper.compile(stanModel, stanHome);
    Command.byPath(compiled)
      .ranIn(Results.getResultFolder())
      .withStandardOutMirroring()
      .withArgs("sample "
          + "thin=" + thin  + " "
          + "num_samples=" + mainOptions.nSamples + " "
          + "num_warmup=" + nWarmUp + " "
          + "random seed=" + mainOptions.random.nextInt()  + " "
          + "data file=" + mainOptions.getRDataFile().getAbsolutePath())
      .call();
  }
}
