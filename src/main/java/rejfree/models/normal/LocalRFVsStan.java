package rejfree.models.normal;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Test;

import com.google.common.base.Stopwatch;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFRunner;
import rejfree.models.normal.NormalChain.NormalChainModel;
import bayonet.coda.EffectiveSize;
import binc.Command;
import blang.variables.RealVariable;
import briefj.BriefFiles;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class LocalRFVsStan implements Runnable
{
  @OptionSet(name = "modelOptions")
  public NormalChainOptions options = new NormalChainOptions();
  
  @OptionSet(name = "rfOptions")
  public RFSamplerOptions rfOptions = new RFSamplerOptions();
  
  @Option
  public int nRepeats = 100;
  
  @Option
  public int nStanIters = 1000;
  
  @Option
  public int nStanWarmUps = 1000;
  
  @Option
  public int essMonitoredVariable = 0;
  
  @Option(required = true)
  public File stanHome = null;
  
  private NormalChain chain;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new LocalRFVsStan());
  }

  @Test
  public void run()
  {
    org.jblas.util.Random.seed(options.random.nextLong());
    chain = new NormalChain(options);
    OutputManager output = new OutputManager();
    output.setOutputFolder(Results.getResultFolder());
    
    File stanProgram = chain.stanProgram(stanHome);
    
    for (int rep = 0; rep < nRepeats; rep++)
    {
      long seed = this.options.random.nextLong();
      DoubleMatrix exactSample = chain.exactSample();
      
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
      int essD = essMonitoredVariable;
      RealVariable aVar = modelSpec.variables.get(essD);
      
      // run stan
      File stanOutput = Results.getFileInResultFolder("currentStanOutput.csv");
      Command stanCommand = createCommand(stanProgram, exactSample, stanOutput);
      Stopwatch watch = Stopwatch.createStarted();
      Command.call(stanCommand);
      long stanRunningTime = watch.elapsed(TimeUnit.MILLISECONDS);
      double stanRunningTimeSec = stanRunningTime / 1000.0;
      Map<String,SummaryStatistics> stanStatistics = stanOutputToSummaryStatistics(stanOutput);
      
      double stanEss = EffectiveSize.effectiveSize(variableSamplesFromStanOutput(stanOutput, stanVarName(essD)));
      Results.getGlobalOutputManager().printWrite("essPerSec", "method", "stan", "value", (stanEss/stanRunningTimeSec));
      
      // run ours for the same time
      LocalRFRunner runner = new LocalRFRunner();
      runner.rfOptions = rfOptions;
      runner.init(modelSpec);
      runner.addMomentRayProcessor();
      runner.addSaveRaysProcessor(Collections.singleton(aVar));
      
      runner.maxSteps = Integer.MAX_VALUE;
      runner.maxTrajectoryLength = Double.POSITIVE_INFINITY;
      runner.maxRunningTimeMilli = stanRunningTime;
      runner.run();
      
      double rfESS = EffectiveSize.effectiveSize(runner.saveRaysProcessor.convertToSample(aVar, 4.0));
      Results.getGlobalOutputManager().printWrite("essPerSec", "method", "RF", "value", (rfESS/stanRunningTimeSec));
      
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable variable = modelSpec.variables.get(d);
        double truth = chain.covarMatrix.get(d, d);
        
        for (boolean isRF : new boolean[]{true,false})
        {
          double estimate = isRF ? runner.momentRayProcessor.getVarianceEstimate(variable) : stanStatistics.get(stanVarName(d)).getVariance();
          double error = Math.abs(truth - estimate);
          output.printWrite("results", 
              "method", (isRF ? "RF" : "STAN"),
              "dim", d, 
              "rep", rep, 
              "seed", seed, 
              "absError", error, 
              "relError", (error/truth), 
              "truth", truth, 
              "estimate", estimate);
        }
      }
    }
    
    output.close();
  }
  
  private String stanVarName(int d)
  {
    return "x." + (d+1);
  }

  private Map<String,SummaryStatistics> stanOutputToSummaryStatistics(File stanOutput)
  {
    Map<String,SummaryStatistics> result = new HashMap<>();
    for (Map<String,String> sample : BriefIO.readLines(stanOutput).indexCSV('#'))
      for (String key : sample.keySet())
        BriefMaps.getOrPut(result, key, new SummaryStatistics()).addValue(Double.parseDouble(sample.get(key)));
    return result;
  }
  
  private List<Double> variableSamplesFromStanOutput(File stanOutput, String varName)
  {
    List<Double> result = new ArrayList<Double>();
    for (Map<String,String> sample : BriefIO.readLines(stanOutput).indexCSV('#'))
      result.add(Double.parseDouble(sample.get(varName)));
    return result;
  }

  private Command createCommand(File stanProgram, DoubleMatrix init, File output)
  {
    File initFile = createInit(init);
    return Command.byPath(stanProgram)
      .ranIn(Results.getResultFolder())
      .withStandardOutMirroring()
      .withArgs(
          "sample " +
            "num_samples=" + nStanIters + " " +
            "num_warmup=" + nStanWarmUps + " " +
          "output " +
            "diagnostic_file=diagnostic.txt" + " " + 
            "file=" + output.getAbsolutePath() + " " +
          "init=" + initFile.getAbsolutePath() + " " +
          "random seed=" + options.random.nextInt());
  }

  private File createInit(DoubleMatrix exactSample)
  {
    File temp = BriefFiles.createTempFile();
    StringBuilder initString = new StringBuilder();
    initString.append("x <- c(");
    for (int d = 0; d < exactSample.length; d++)
      initString.append("" + exactSample.get(d) + (d == exactSample.length - 1 ? "" : ","));
    initString.append(")");
    BriefIO.write(temp, initString);
    return temp;
  }

}
