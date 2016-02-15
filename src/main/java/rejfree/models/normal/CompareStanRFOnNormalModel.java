package rejfree.models.normal;

import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Test;

import rejfree.RFSamplerOptions;
import rejfree.StanUtils;
import rejfree.StanUtils.StanExecution;
import rejfree.local.LocalRFRunner;
import rejfree.local.TrajectoryRay;
import rejfree.models.normal.NormalChain.NormalChainModel;
import rejfree.processors.MomentRayProcessor;
import rejfree.processors.SaveRaysProcessor;
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
      
      Map<String, List<Double>> variableSamplesFromStanOutput = needStan ? stanExec.parsedStanOutput() : null;
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
            List<Double> stanSamples = isRF 
                ? null //runner.saveRaysProcessor.convertToSample(variable, 0.01)
                : variableSamplesFromStanOutput.get(stanVarName(d));
            if (recordPartialSums)
            {
              if (isRF)
                recordPartialSums(runner.saveRaysProcessor, methodName, d, truth, variable, time, rep);
              else
                recordPartialSums(stanSamples, methodName, d, truth, time, rep);
            }
            
            // commented because seems very sensitive to sampling rate
//            double ess = EffectiveSize.effectiveSize(samples);
//            Results.getGlobalOutputManager().printWrite("essPerSec", 
//                "method", methodName, 
//                "dim", d,
//                "rep", rep,
//                "value", (1000.0*ess/runner.options.maxRunningTimeMilli));
            
            if (fractionOfTrajectoryToPlot != 0.0)
              plotTrajectory(isRF ? runner.saveRaysProcessor.convertToSample(variable, 1) : stanSamples, output, isRF ? "RF" : "STAN", d);
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
  



  private void recordPartialSums(SaveRaysProcessor saveRaysProcessor,
      String methodName, int d, double truth, RealVariable v, double fullTime, int rep)
  {
    List<TrajectoryRay> rays = saveRaysProcessor.samples.get(v);
    double sumSq = 0.0;
    double sum = 0.0;
    final int interval = rays.size() / 100;
    for (int i = 0; i < rays.size(); i++)  
    {
      TrajectoryRay 
        curRay = rays.get(i),
        nxtRay = i == rays.size() - 1 ? null : rays.get(i+1);
      double endTimeForRay = (nxtRay == null ? saveRaysProcessor.totalLength : nxtRay.t);
      double rayLen = endTimeForRay - curRay.t;
      sumSq += MomentRayProcessor.indefIntegralForVar (curRay.position_t, curRay.velocity_t, rayLen);
      sum +=   MomentRayProcessor.indefIntegralForMean(curRay.position_t, curRay.velocity_t, rayLen);

      if (i > 0 && i % interval == 0)
      {
        final double curTime = endTimeForRay;
        double meanEstimate = sum / curTime;
        final double varEstimate =  sumSq / curTime - (meanEstimate * meanEstimate);
        double error = Math.abs(truth - varEstimate);
        Results.getGlobalOutputManager().printWrite("partialSums",
            "method", methodName,
            "dim", d,
            "percent", (i/interval),
            "rep", rep, 
            "wallClockTime", (i/interval)*fullTime/100,
            "absError", error, 
            "relError", (error/truth), 
            "truth", truth, 
            "estimate", varEstimate);
      }
    }
  }

  private void plotTrajectory(List<Double> samples, OutputManager output, String method, int d)
  {
    double N = (int) (samples.size() * fractionOfTrajectoryToPlot);
    
    for (double i = 0; i < nTrajectoryPoints; i++)
      output.write("trajectories", "step", (double) (100 * i * fractionOfTrajectoryToPlot / nTrajectoryPoints), "method", method, "d", d, "value", samples.get((int) (i * N / nTrajectoryPoints)));
  }
  


  private void recordPartialSums(List<Double> samples, String methodName, int d, double truth, double fullTime, int rep)
  {
    SummaryStatistics stat = new SummaryStatistics();
    int counter = 0;
    final int interval = samples.size() / 100;
    for (double cur : samples)
    {
      stat.addValue(cur);
      if (counter > 0 && counter % interval == 0)
      {
        final double varEstimate =  stat.getVariance();
        double error = Math.abs(truth - varEstimate);
        Results.getGlobalOutputManager().printWrite("partialSums",
            "method", methodName,
            "dim", d,
            "percent", (counter/interval),
            "rep", rep, 
            "wallClockTime", (counter/interval)*fullTime/100,
            "absError", error, 
            "relError", (error/truth), 
            "truth", truth, 
            "estimate", varEstimate);
      }
      counter++;
    }
  }

  public static final String VAR_NAME = "x";
  
  public static String stanVarName(int d)
  {
    return VAR_NAME + "." + (d+1);
  }
}
