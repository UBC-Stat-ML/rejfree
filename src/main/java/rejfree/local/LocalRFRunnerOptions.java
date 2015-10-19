package rejfree.local;

import java.util.Random;

import rejfree.global.GlobalRFSampler.RFSamplerOptions;
import briefj.opt.Option;
import briefj.opt.OptionSet;



public class LocalRFRunnerOptions
{
  @Option
  public long maxRunningTimeMilli = Long.MAX_VALUE;
  
  @Option
  public int maxSteps = 1000;
  
  @Option
  public double maxTrajectoryLength = Double.POSITIVE_INFINITY;
  
  @Option
  public Random samplingRandom = new Random(1);
  
  @OptionSet(name = "rfOptions")
  public RFSamplerOptions rfOptions = new RFSamplerOptions();
}