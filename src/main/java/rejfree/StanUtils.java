package rejfree;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;

import com.google.common.base.Stopwatch;

import bayonet.bugs.StanWrapper;
import binc.Command;
import briefj.BriefFiles;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.opt.Option;
import briefj.run.Results;



public class StanUtils
{
  public static class StanOptions
  {
    @Option
    public int nStanIters = 1000;
    
    @Option
    public int nStanWarmUps = 1000;
    
    @Option
    public File stanHome = null;
    
    @Option 
    public boolean useNuts = true;
    
    @Option
    public String outputName = "stan-output.csv";
    
    @Option
    public Random rand = new Random(1);

    @Option
    public boolean useDiagMetric = true;

    @Option
    public double stepSize = 1.0;
    
    @Option
    public boolean verbose = false;
    
    @Option
    public double stepSizeJitter = 0.0;
    
    @Option
    public double intTime = 2.0 * Math.PI;
  }
  
  public static class StanExecution
  {
    public final StanOptions options;
    public final String model;
    private long elapsed;
    private File output;
    private boolean ran = false;
    private StringBuilder initString = new StringBuilder();
    private File dataFile;
      
    public StanExecution(String model, StanOptions options)
    {
      this(model, options, null);
    }
    
    public StanExecution(String model, StanOptions options, File dataFile)
    {
      this.dataFile = dataFile;
      this.options = options;
      this.model = model;
      output = Results.getFileInResultFolder(options.outputName);
    }
    
    public void addInit(String varName, DoubleMatrix values)
    {
      initString.append(varName + " <- c(");
      for (int d = 0; d < values.length; d++)
        initString.append("" + values.get(d) + (d == values.length - 1 ? "" : ","));
      initString.append(")\n");
    }
    
    private boolean hasDataFile()
    {
      return dataFile != null;
    }

    private File stanProgram()
    {
      return StanWrapper.compile(model, options.stanHome);
    }
    
    public void run()
    {
      Command stanCommand = createCommand();
      Stopwatch watch = Stopwatch.createStarted();
      Command.call(stanCommand);
      elapsed = watch.elapsed(TimeUnit.MILLISECONDS);
      ran = true;
    }
    
    public long getRunningTimeMilli()
    {
      checkRan();
      return elapsed;
    }
    
    private void checkRan()
    {
      if (!ran)
        throw new RuntimeException("Need to call run() first");
    }
    
    private Command createCommand()
    {
      boolean hasInit = !initString.toString().isEmpty();
      File initFile = hasInit ? createInit() : null;
      File stanProgram = stanProgram();
      if (options.verbose)
        System.out.println("Running stan exec cached at:" + stanProgram.getAbsolutePath());
      Command result =  Command.byPath(stanProgram)
        .ranIn(Results.getResultFolder())
        .withArgs(
            "sample " +
              "num_samples=" + options.nStanIters + " " +
              "num_warmup=" + options.nStanWarmUps + " " +
              "adapt " +
                "engaged=" + (options.nStanWarmUps == 0 ? "0" : "1") + " " +
              "algorithm=hmc " +
                "engine=" + (options.useNuts ? "nuts" : "static") + " " +
                  (options.useNuts ? "" : "int_time=" + options.intTime) + " " +
                "metric=" + (options.useDiagMetric  ? "diag_e" : "unit_e") + " " +
                "stepsize=" + options.stepSize + " " + 
                "stepsize_jitter=" + options.stepSizeJitter + " " + 
            (hasDataFile() ? "data file=" + dataFile.getAbsolutePath() + " " : "") + 
            "output " +
              "file=" + output.getAbsolutePath() + " " +
            (hasInit ? "init=" + initFile.getAbsolutePath() : "") + " " +
            "random seed=" + options.rand.nextInt());
      if (options.verbose)
        result = result.withStandardOutMirroring();
      return result;
    }
    
    public Map<String,SummaryStatistics> stanOutputToSummaryStatistics()
    {
      checkRan();
      Map<String,SummaryStatistics> result = new HashMap<>();
      for (Map<String,String> sample : BriefIO.readLines(output).indexCSV('#'))
        for (String key : sample.keySet())
          BriefMaps.getOrPut(result, key, new SummaryStatistics()).addValue(Double.parseDouble(sample.get(key)));
      return result;
    }
    
    public Map<String,List<Double>> variableSamplesFromStanOutput()
    {
      Map<String,List<Double>> result = new LinkedHashMap<String, List<Double>>();
      for (Map<String,String> sample : BriefIO.readLines(output).indexCSV('#'))
        for (String varName : sample.keySet())
          BriefMaps.getOrPutList(result, varName).add(Double.parseDouble(sample.get(varName)));
      return result;
    }
    
    private File createInit()
    {
      File temp = BriefFiles.createTempFile();
      BriefIO.write(temp, initString.toString());
      return temp;
    }
  }
}
