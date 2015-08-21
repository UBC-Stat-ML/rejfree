package rejfree.models.normal;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Test;
import org.mvel2.templates.TemplateRuntime;

import com.google.common.base.Stopwatch;

import rejfree.global.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFRunner;
import bayonet.bugs.StanWrapper;
import bayonet.coda.EffectiveSize;
import binc.Command;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.BriefFiles;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class CompareStanRFOnRN implements Runnable
{
  @Option
  public int nVars = 100;
  
  @Option
  public Random rand = new Random(1);
  
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
  
  @Option
  public boolean useLocal = true;
  
  @Option(required = true)
  public File stanHome = null;
  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new CompareStanRFOnRN());
  }
  
  private String stanModel()
  {
    String template = BriefIO.resourceToString("/rejfree/stanRNTemplate.txt");
    return (String) TemplateRuntime.eval(template, this);
  }
  
  public File stanProgram(File stanHome)
  {
    return StanWrapper.compile(stanModel(), stanHome);
  }

  @Test
  public void run()
  {
    org.jblas.util.Random.seed(rand.nextLong());
    
    OutputManager output = new OutputManager();
    output.setOutputFolder(Results.getResultFolder());
    
    File stanProgram = stanProgram(stanHome);
    
    
    
    for (int rep = 0; rep < nRepeats; rep++)
    {
      long seed = rand.nextLong();
      DoubleMatrix exactSample = new DoubleMatrix(nVars);
      
      if (!useLocal)
        throw new RuntimeException();
      
      NormalRNModel2 modelSpec = 
          //useLocal ? 
          new NormalRNModel2(exactSample.data) 
      ;
      //: new NormalRNModel(exactSample.data);
      RealVariable aVar = modelSpec.variables.get(0);
      
      int essD = essMonitoredVariable;
      
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
      runner.options.rfOptions = rfOptions;
      runner.init(modelSpec);
      runner.addMomentRayProcessor();
      runner.addSaveRaysProcessor(Collections.singleton(aVar));
      
      runner.options.maxSteps = Integer.MAX_VALUE;
      runner.options.maxTrajectoryLength = Double.POSITIVE_INFINITY;
      runner.options.maxRunningTimeMilli = stanRunningTime;
      runner.run();
      
      double rfESS = EffectiveSize.effectiveSize(runner.saveRaysProcessor.convertToSample(aVar, 4.0));
      Results.getGlobalOutputManager().printWrite("essPerSec", "method", "RF", "value", (rfESS/stanRunningTimeSec));
      
      double delta = 1.0 / ((double) nVars);
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable variable = modelSpec.variables.get(d);
        
        for (boolean isRF : new boolean[]{true,false})
        {
          double estimate = isRF ? runner.momentRayProcessor.getVarianceEstimate(variable) : stanStatistics.get(stanVarName(d)).getVariance();
          double stdDev = ((double) (d+1)) * delta;
          double truth = stdDev * stdDev;
          
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
  

  
  public class NormalRNModel2 
  {
    public List<RealVariable> variables = new ArrayList<>();
    
    @DefineFactor
    public List<UnivariateNormalFactor> localFactors;
    
    public NormalRNModel2(double [] init)
    {
      this.localFactors = init(init) ;
    }
    
    private List<UnivariateNormalFactor> init(double [] init)
    {
      double delta = 1.0 / ((double) nVars);
      List<UnivariateNormalFactor> result = new ArrayList<>();
      for (int i = 0; i < nVars; i++)
      {
        RealVariable cur = RealVariable.real(init[i]);
        variables.add(cur);
        double stdDev = ((double) (i+1)) * delta;
        double variance = stdDev * stdDev;
        UnivariateNormalFactor curFac = new UnivariateNormalFactor(cur, variance);
        result.add(curFac);
      }
      return result;
    }
  }
  
  public class NormalRNModel
  {
    public List<RealVariable> variables = new ArrayList<>();
    
    @DefineFactor
    public RNNormalFactor localFactor;
    
    public NormalRNModel(double [] init)
    {
      this.localFactor = init(init) ;
    }
    
    private RNNormalFactor init(double [] init)
    {
      for (int i = 0; i < nVars; i++)
      {
        RealVariable cur = RealVariable.real(init[i]);
        variables.add(cur);
      }
      return new RNNormalFactor(variables);
    }

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
          "random seed=" + rand.nextInt());
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
