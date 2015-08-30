package rejfree.models.normal;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.junit.Test;
import org.mvel2.templates.TemplateRuntime;

import rejfree.StanUtils;
import rejfree.global.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFRunner;
import bayonet.coda.EffectiveSize;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.BriefIO;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class CompareStanRFOnRN implements Runnable
{
  @Option
  public int nVars = 100;
  
  @OptionSet(name = "rfOptions")
  public RFSamplerOptions rfOptions = new RFSamplerOptions();
  
  @OptionSet(name = "stan")
  public StanUtils.StanOptions stanOptions = new StanUtils.StanOptions(); 
  
  @Option
  public int nRepeats = 100;
  
  @Option(gloss = "Which variables to monitor (1 for monitoring all of them)")
  public int variableMonitorInterval = 1;
  
  @Option
  public boolean useLocal = true;
  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new CompareStanRFOnRN());
  }
  
  private String stanModel()
  {
    String template = BriefIO.resourceToString("/rejfree/stanRNTemplate.txt");
    return (String) TemplateRuntime.eval(template, this);
  }
  
  @Test
  public void run()
  {
    OutputManager output = Results.getGlobalOutputManager();
    
    for (int rep = 0; rep < nRepeats; rep++)
    {
      DoubleMatrix exactSample = new DoubleMatrix(nVars);
      
      BaseModel modelSpec = 
          useLocal ? 
          new NormalRNModelLocal(exactSample.data) :
          new NormalRNModelGlobal(exactSample.data);
      
      // run stan
      StanUtils.StanExecution stanExec = new StanUtils.StanExecution(stanModel(), stanOptions);
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
      
      double delta = 1.0 / ((double) nVars);
      for (int d = 0; d < modelSpec.variables.size(); d++)
      {
        RealVariable variable = modelSpec.variables.get(d);
        
        
        for (boolean isRF : new boolean[]{true,false})
        {
          double estimate = isRF ? runner.momentRayProcessor.getVarianceEstimate(variable) : stanStatistics.get(stanVarName(d)).getVariance();
          double stdDev = ((double) (d+1)) * delta;
          double truth = stdDev * stdDev;
          final String methodName = (isRF ? "RF" : "STAN");
          double error = Math.abs(truth - estimate);
          output.printWrite("results", 
              "method", methodName,
              "dim", d, 
              "rep", rep, 
              "absError", error, 
              "relError", (error/truth), 
              "truth", truth, 
              "estimate", estimate);
          
          if (d % variableMonitorInterval == 0)
          {
            double ess = EffectiveSize.effectiveSize( isRF 
                ? runner.saveRaysProcessor.convertToSample(variable, 4.0)
                : variableSamplesFromStanOutput.get(stanVarName(d)));
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
  
  public class BaseModel
  {
    public List<RealVariable> variables = new ArrayList<>();
  }
  
  public class NormalRNModelLocal extends BaseModel
  {
    @DefineFactor
    public List<UnivariateNormalFactor> localFactors;
    
    public NormalRNModelLocal(double [] init)
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
  
  public class NormalRNModelGlobal extends BaseModel
  {
    
    @DefineFactor
    public RNNormalFactor localFactor;
    
    public NormalRNModelGlobal(double [] init)
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
  
  public static final String VAR_NAME = "x";
  
  private String stanVarName(int d)
  {
    return VAR_NAME + "." + (d+1);
  }
}
