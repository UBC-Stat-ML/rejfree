package rejfree.models.normal;

import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Test;

import rejfree.StanUtils;
import rejfree.StanUtils.StanExecution;
import rejfree.global.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFRunner;
import rejfree.models.normal.NormalChain.NormalChainModel;
import bayonet.coda.EffectiveSize;
import bayonet.rplot.RUtils;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class CompareStanRFOnNormalModel implements Runnable
{
  @OptionSet(name = "modelOptions")
  public NormalChainOptions options = new NormalChainOptions();
  
  @OptionSet(name = "rfOptions")
  public RFSamplerOptions rfOptions = new RFSamplerOptions();
  
  @OptionSet(name = "stan")
  public StanUtils.StanOptions stanOptions = new StanUtils.StanOptions();
  
  @Option(gloss = "If equal to zero, run stan first to compare against it.")
  public long timeMilli = 0;
  
  @Option 
  public boolean recordPartialSums = false;
  
  @Option
  public int nRepeats = 100;
  
  @Option(gloss = "Which variables to monitor (1 for monitoring all of them)")
  public int variableMonitorInterval = 10;
  
  @Option 
  public double fractionOfTrajectoryToPlot = 0.0;
  
  @Option
  public boolean saveTrajectories = false;
  
  @Option 
  public int nTrajectoryPoints = 10000;
  
  private NormalChain chain;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new CompareStanRFOnNormalModel());
  }
  
  @Test
  public void run()
  {
    chain = new NormalChain(options);
    OutputManager output = Results.getGlobalOutputManager();
    
    boolean needStan = (timeMilli == 0);
    boolean[] isRFs = needStan ? new boolean[]{true,false} : new boolean[]{true};
    
    for (int rep = 0; rep < nRepeats; rep++)
    {
      DoubleMatrix exactSample = chain.exactSample();
      
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
      
      long time = timeMilli;
      Map<String,SummaryStatistics> stanStatistics = null;
      StanExecution stanExec = null;
      if (needStan)
      {
        System.out.println("Running Stan to determine running time cap for RF");
        stanExec = new StanExecution(chain.stanModel(), stanOptions);
        stanExec.addInit(VAR_NAME, exactSample);
        stanExec.run();
        stanStatistics = stanExec.stanOutputToSummaryStatistics();
        time = (long) ((double) stanExec.getRunningTimeMilli() * 0.99);
      }
      else
        System.out.println("Running RF only (skipping Stan) since a running time has been specified. \nIf you want to run Stan leave timeMilli to zero.");
      
      // run ours for the same time
      LocalRFRunner runner = new LocalRFRunner();
      runner.options.rfOptions = rfOptions;
      runner.init(modelSpec);
      runner.addMomentRayProcessor();
      runner.addSaveRaysProcessor(ComparisonUtils.subsample(modelSpec.variables, variableMonitorInterval));
      
      runner.options.maxSteps = Integer.MAX_VALUE;
      runner.options.maxTrajectoryLength = Double.POSITIVE_INFINITY;
      runner.options.maxRunningTimeMilli = time;
      runner.run();
      
      Map<String, List<Double>> variableSamplesFromStanOutput = needStan ? stanExec.variableSamplesFromStanOutput() : null;
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable variable = modelSpec.variables.get(d);
        double truth = chain.covarMatrix.get(d, d);
        
        for (boolean isRF : isRFs)
        {
          double estimate = isRF ? runner.momentRayProcessor.getVarianceEstimate(variable) : stanStatistics.get(stanVarName(d)).getVariance();
          double error = Math.abs(truth - estimate);
          final String methodName = (isRF ? "RF" : "STAN");
          output.printWrite("varianceEstimates", 
              "method", methodName,
              "dim", d, 
              "rep", rep, 
              "absError", error, 
              "relError", (error/truth), 
              "truth", truth, 
              "estimate", estimate);
          
          if (d % variableMonitorInterval == 0)
          {
            List<Double> samples = isRF 
                ? runner.saveRaysProcessor.convertToSample(variable, 4.0)
                : variableSamplesFromStanOutput.get(stanVarName(d));
            if (recordPartialSums)
              recordPartialSums(samples, methodName, d, truth);
            double ess = EffectiveSize.effectiveSize(samples);
            Results.getGlobalOutputManager().printWrite("essPerSec", 
                "method", methodName, 
                "dim", d,
                "rep", rep,
                "value", (1000.0*ess/runner.options.maxRunningTimeMilli));
            
            if (fractionOfTrajectoryToPlot != 0.0)
              plotTrajectory(samples, output, isRF ? "RF" : "STAN", d);
          }
        }
      }
      
      
      if (saveTrajectories)
        runner.saveRaysProcessor.toCSV(Results.getFileInResultFolder("full-rf-trajectory-" + rep + ".csv"));
    }
    
    output.close();
    
    if (fractionOfTrajectoryToPlot != 0.0)
      RUtils.callGeneratedRScript("/rejfree/plotTrajectory2.txt", Pair.of(Results.getFileInResultFolder("trajectories.csv"), Results.getFileInResultFolder("trajectories.pdf")));
  }
  
  private void plotTrajectory(List<Double> samples, OutputManager output, String method, int d)
  {
    double N = (int) (samples.size() * fractionOfTrajectoryToPlot);
    
    for (double i = 0; i < nTrajectoryPoints; i++)
      output.write("trajectories", "step", (double) (100 * i * fractionOfTrajectoryToPlot / nTrajectoryPoints), "method", method, "d", d, "value", samples.get((int) (i * N / nTrajectoryPoints)));
  }

  private void recordPartialSums(List<Double> samples, String methodName, int dim, double truth)
  {
    SummaryStatistics stat = new SummaryStatistics();
    int counter = 0;
    final int interval = samples.size() / 100;
    int percent = 0;
    for (double cur : samples)
    {
      stat.addValue(cur);
      if (counter > 0 && counter % interval == 0)
        Results.getGlobalOutputManager().printWrite("partialSums",
            "method", methodName,
            "dim", dim,
            "percent", percent++,
            "value", (stat.getVariance() - truth) /truth);
      counter++;
    }
  }

  public static final String VAR_NAME = "x";
  
  private String stanVarName(int d)
  {
    return VAR_NAME + "." + (d+1);
  }
}
