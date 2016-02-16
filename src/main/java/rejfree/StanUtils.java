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
    public boolean saveWarmUp = false;
    
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
    public final File output;
    private boolean ran = false;
    private final StringBuilder initString = new StringBuilder();
    private final StringBuilder dataString = new StringBuilder();
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
      add(varName, values, initString);
    }
    
    public void addInit(String varName, Number value)
    {
      add(varName, value, initString);
    }
    
    public void addData(String varName, DoubleMatrix values)
    {
      add(varName, values, dataString);
    }
    
    public void addData(String varName, Number value)
    {
      add(varName, value, dataString);
    }
    
    private static void add(String varName, DoubleMatrix values, StringBuilder toAddTo)
    {
      toAddTo.append(varName + " <- c(");
      for (int d = 0; d < values.length; d++)
        toAddTo.append("" + values.get(d) + (d == values.length - 1 ? "" : ","));
      toAddTo.append(")\n");
    }
    
    private static void add(String varName, Number value, StringBuilder toAddTo)
    {
      toAddTo.append(varName + " <- " + value + "\n");
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
      File dataFile = createData();
      Command result =  Command.byPath(stanProgram)
        .ranIn(Results.getResultFolder())
        .withArgs(
            "sample " +
              "num_samples=" + options.nStanIters + " " +
              "num_warmup=" + options.nStanWarmUps + " " +
              "save_warmup=" + (options.saveWarmUp ? "1" : "0") + " " + // could not get this to work
              "adapt " +
                "engaged=" + (options.nStanWarmUps == 0 ? "0" : "1") + " " +
              "algorithm=hmc " +
                "engine=" + (options.useNuts ? "nuts" : "static") + " " +
                  (options.useNuts ? "" : "int_time=" + options.intTime) + " " +
                "metric=" + (options.useDiagMetric  ? "diag_e" : "unit_e") + " " +
                "stepsize=" + options.stepSize + " " + 
                "stepsize_jitter=" + options.stepSizeJitter + " " + 
            (dataFile != null ? "data file=" + dataFile.getAbsolutePath() + " " : "") + 
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
    
    /**
     * Contains the value of the variable at each iteration and moreover statistics such as acceptance rate, log probab, etc.
     */
    public Map<String,List<Double>> parsedStanOutput()
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
    
    private File createData()
    {
      if (dataFile == null && dataString.length() == 0)
        return null;
      if (dataFile == null)
        dataFile = BriefFiles.createTempFile();
      else
        dataString.append(BriefIO.fileToString(dataFile));
      BriefIO.write(dataFile, dataString.toString());
      return dataFile;
    }
  }
}
