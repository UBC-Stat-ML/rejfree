package rejfree.local;

import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import com.google.common.base.Stopwatch;

import rejfree.global.GlobalRFSampler.RFSamplerOptions;
import rejfree.processors.MomentRayProcessor;
import rejfree.processors.SaveRaysProcessor;
import rejfree.processors.SaveSamplesProcessor;
import blang.ProbabilityModel;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Results;



public class LocalRFRunner
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
  
  public ProbabilityModel model;
  public LocalRFSampler sampler;
  
  private OutputManager output = Results.getGlobalOutputManager();
  
  public MomentRayProcessor momentRayProcessor = null;
  public SaveRaysProcessor saveRaysProcessor = null;
  public SaveSamplesProcessor saveSamplesProcessor = null;
  
  public void init(Object modelSpec)
  {
    if (isInit())
      throw new RuntimeException("Model already initialized.");
    model = new ProbabilityModel(modelSpec);
    sampler = new LocalRFSampler(model, rfOptions);
  }
  
  private boolean isInit()
  {
    return model != null;
  }
  
  private void checkInit()
  {
    if (!isInit())
      throw new RuntimeException("Model should first be initialized.");
  }
  
  public void addMomentRayProcessor()
  {
    checkInit();
    if (momentRayProcessor != null)
      throw new RuntimeException("Alreayd added");
    momentRayProcessor = new MomentRayProcessor();
    sampler.addRayProcessor(momentRayProcessor);
  }
  
  public void addSaveRaysProcessor(Collection<RealVariable> variables)
  {
    checkInit();
    if (saveRaysProcessor != null)
      throw new RuntimeException("Alreayd added");
    saveRaysProcessor = new SaveRaysProcessor(variables);
    sampler.addRayProcessor(saveRaysProcessor);
  }
  
  public void addSaveAllRaysProcessor()
  {
    checkInit();
    List<RealVariable> all = model.getLatentVariables(RealVariable.class);
    addSaveRaysProcessor(all);
  }
  
  public void addSaveSamplesProcessor(List<RealVariable> variables)
  {
    checkInit();
    if (saveSamplesProcessor != null)
      throw new RuntimeException("Alreayd added");
    saveSamplesProcessor = new SaveSamplesProcessor(variables);
    sampler.addPointProcessor(saveSamplesProcessor);
  }

  public void run()
  {
    checkInit();
    
    Stopwatch watch = Stopwatch.createStarted();
    sampler.iterate(samplingRandom, maxSteps, maxTrajectoryLength, maxRunningTimeMilli);
    long elapsed = watch.elapsed(TimeUnit.MILLISECONDS);
    
    output.printWrite("general-sampler-diagnostic", 
        "wallClockTimeMilli", elapsed, 
        "nCollisions", sampler.getNCollisions(),
        "nCollidedVariables", sampler.getNCollidedVariables(),
        "nRefreshments", sampler.getNRefreshments(),
        "nRefreshedVariables", sampler.getNRefreshedVariables());
    
    output.flush();
  }


}
