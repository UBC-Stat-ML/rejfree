package rejfree.models.expfam;

import java.io.File;
import java.util.Random;

import rejfree.StanUtils;
import rejfree.StanUtils.StanExecution;
import rejfree.local.LocalRFRunner;
import rejfree.local.LocalRFRunnerOptions;
import rejfree.models.expfam.MRF.ModelSpec;
import briefj.BriefIO;
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
    
    
    if (runStan)
    {
      System.out.println("Starting Stan");
      ModelSpec modelSpec = mrf.newModelSpecFromGenerativeProcess(new Random(generateRandom));
      File stanData = Results.getFileInResultFolder("stanDataString.txt");
      BriefIO.write(stanData, modelSpec.getDataString());
      StanExecution stanExec = new StanExecution(mrf.stanModel(), stanOptions, stanData);
      stanExec.addInit(mrf.stanLatentVariableName(), modelSpec.getLatentAsDoubleMatrix());
      stanExec.run();
    }
    
    if (runRF)
    {
      System.out.println("Starting RF");
      ModelSpec modelSpec = mrf.newModelSpecFromGenerativeProcess(new Random(generateRandom));
      LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);
      rfRunner.init(modelSpec);
      rfRunner.run();
    }
  }

}
