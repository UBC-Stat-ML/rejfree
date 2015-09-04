package rejfree.models.normal;

import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Test;

import rejfree.StanUtils;
import rejfree.StanUtils.StanExecution;
import rejfree.global.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFRunner;
import rejfree.models.normal.NormalChain.NormalChainModel;
import bayonet.coda.EffectiveSize;
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
  
  @Option
  public int nRepeats = 100;
  
  @Option(gloss = "Which variables to monitor (1 for monitoring all of them)")
  public int variableMonitorInterval = 10;
  
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
    
    for (int rep = 0; rep < nRepeats; rep++)
    {
      DoubleMatrix exactSample = chain.exactSample();
      
      NormalChainModel modelSpec = chain.new NormalChainModel(exactSample.data);
      
      StanExecution stanExec = new StanExecution(chain.stanModel(), stanOptions);
      stanExec.addInit(VAR_NAME, exactSample);
      stanExec.run();
      Map<String,SummaryStatistics> stanStatistics = stanExec.stanOutputToSummaryStatistics();
      
      // run ours for the same time
      LocalRFRunner runner = new LocalRFRunner();
      runner.options.rfOptions = rfOptions;
      runner.init(modelSpec);
      runner.addMomentRayProcessor();
      runner.addSaveRaysProcessor(ComparisonUtils.subsample(modelSpec.variables, variableMonitorInterval));
      
      runner.options.maxSteps = Integer.MAX_VALUE;
      runner.options.maxTrajectoryLength = Double.POSITIVE_INFINITY;
      runner.options.maxRunningTimeMilli = (long) ((double) stanExec.getRunningTimeMilli() * 0.99);
      runner.run();
      
      Map<String, List<Double>> variableSamplesFromStanOutput = stanExec.variableSamplesFromStanOutput();
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable variable = modelSpec.variables.get(d);
        double truth = chain.covarMatrix.get(d, d);
        
        for (boolean isRF : new boolean[]{true,false})
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
            recordPartialSums(samples, methodName, d, truth);
            double ess = EffectiveSize.effectiveSize(samples);
            Results.getGlobalOutputManager().printWrite("essPerSec", 
                "method", methodName, 
                "dim", d,
                "rep", rep,
                "value", (1000.0*ess/runner.options.maxRunningTimeMilli));
          }
        }
      }
    }
    
    output.close();
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
