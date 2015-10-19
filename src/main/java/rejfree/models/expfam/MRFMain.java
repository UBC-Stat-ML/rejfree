package rejfree.models.expfam;

import java.io.File;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import rejfree.StanUtils;
import rejfree.StanUtils.StanExecution;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.models.expfam.MRF.ModelSpec;
import rejfree.processors.MomentRayProcessor;
import briefj.BriefIO;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class MRFMain implements Runnable
{
  @OptionSet(name = "rfOptions")
  public LocalRFRunnerOptions localRFRunnerOption = new LocalRFRunnerOptions();
  
  @OptionSet(name = "stan")
  public StanUtils.StanOptions stanOptions = new StanUtils.StanOptions(); 
  
  @OptionSet(name = "model")
  public MRFOptions modelOptions = new MRFOptions();
  
  @Option
  public boolean runRF = true;
  
  @Option
  public boolean runStan = true;
  
  @Option
  public long generateRandom = 31;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new MRFMain());
  }

  @Override
  public void run()
  {
    MRF mrf = new MRF(modelOptions);
    OutputManager output = Results.getGlobalOutputManager();
    
    long stanRuntimeMillis = -1;
    if (runStan)
    {
      System.out.println("Starting Stan");
      ModelSpec modelSpec = mrf.newModelSpecFromGenerativeProcess(new Random(generateRandom));
      File stanData = Results.getFileInResultFolder("stanDataString.txt");
      BriefIO.write(stanData, modelSpec.getDataString());
      StanExecution stanExec = new StanExecution(mrf.stanModel(), stanOptions, stanData);
      stanExec.addInit(mrf.stanLatentVariableName(), modelSpec.getLatentAsDoubleMatrix());
      stanExec.run();
      stanRuntimeMillis = stanExec.getRunningTimeMilli();
      Map<String, SummaryStatistics> stanOutputToSummaryStatistics = stanExec.stanOutputToSummaryStatistics();
      
      for (int var = 0; var < mrf.nLatentVariables(); var++)
      {
        output.printWrite("moments", "method", "Stan", "var", var, "stat", "mean", "value", stanOutputToSummaryStatistics.get(mrf.stanLatentCoordinate(var)).getMean());
        output.printWrite("moments", "method", "Stan", "var", var, "stat", "var",  "value", stanOutputToSummaryStatistics.get(mrf.stanLatentCoordinate(var)).getVariance());
      }
    }
    
    if (runRF)
    {
      if (runStan)
      {
        localRFRunnerOption.maxRunningTimeMilli = stanRuntimeMillis;
        localRFRunnerOption.maxSteps = Integer.MAX_VALUE;
        localRFRunnerOption.maxTrajectoryLength = Double.POSITIVE_INFINITY;
      }
      System.out.println("Starting RF");
      ModelSpec modelSpec = mrf.newModelSpecFromGenerativeProcess(new Random(generateRandom));
      LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);
      rfRunner.init(modelSpec);
      rfRunner.addMomentRayProcessor();
      rfRunner.run();
      MomentRayProcessor rfMomentProcessor = rfRunner.momentRayProcessor;
      for (int var = 0; var < mrf.nLatentVariables(); var++)
      {
        output.printWrite("moments", "method", "RF", "var", var, "stat", "mean", "value", rfMomentProcessor.getMeanEstimate(modelSpec.latentVariables.get(var)));
        output.printWrite("moments", "method", "RF", "var", var, "stat", "var",  "value", rfMomentProcessor.getVarianceEstimate(modelSpec.latentVariables.get(var)));
      }
    }
    
    output.close();
  }

}
