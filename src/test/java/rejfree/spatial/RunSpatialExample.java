package rejfree.spatial;

import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class RunSpatialExample implements Runnable
{
  @OptionSet(name = "main")
  public SpatialMainOptions mainOptions = new SpatialMainOptions();
  
  @OptionSet(name = "rf")
  public SpatialBlang rfMain = new SpatialBlang(mainOptions);
  
  @OptionSet(name = "hmc")
  public SpatialStan stanMain = new SpatialStan(mainOptions);
  
  @Option
  public SamplingMethod samplingMethod = SamplingMethod.LOCAL_RF;
  
  public static enum SamplingMethod
  {
    STAN     { void execute(RunSpatialExample instance) { instance.stanMain.run(); }},
    LOCAL_RF { void execute(RunSpatialExample instance) { instance.rfMain.run(); }};
    abstract void execute(RunSpatialExample instance);
  }

  @Override
  public void run()
  {
    // run
    samplingMethod.execute(this);
    
    // post-process
    new PostprocessIntersectionSamples(Results.getFileInResultFolder(SAMPLES_FILE_NAME)).run();
  }
  
  // NB: do not change this value (we rely this to be the same as stan's default)
  public static String SAMPLES_FILE_NAME = "output.csv";
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new RunSpatialExample());
  }
}
