package rejfree.local;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import rejfree.local.CompareESSLocalGlobal.ModelSpec;
import bayonet.rplot.PlotLine;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class CheckMCErrorScaling implements Runnable
{
  @OptionSet(name = "rf")
  public LocalRFRunnerOptions options = new LocalRFRunnerOptions();
  
  @Option
  public int nRepeats = 100;

  @Option
  public boolean useSparse = false;
  
  @Option
  public boolean exact = false;

  @Option
  public Random initRandom = new Random(1);
  
  @Option
  public int dim = 1;
  
  private OutputManager output = Results.getGlobalOutputManager();

  @Option
  public int minN = 100;
  
  @Option
  public int maxN = 10000;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new CheckMCErrorScaling());
  }
  
  public double mcError(int N)
  {
    options.maxSteps = Integer.MAX_VALUE;
    options.maxTrajectoryLength = N;
    SummaryStatistics mcErrorStatistics = new SummaryStatistics();
    ModelSpec spec = new CompareESSLocalGlobal.ModelSpec(dim, useSparse);
    spec.initFromStatio(initRandom);
    for (int i = 0; i < nRepeats; i++)
    {
      double estimate = Double.NaN;
      if (exact)
      {
        SummaryStatistics exact = new SummaryStatistics();
        for (int j = 0; j < N; j++)
          exact.addValue(initRandom.nextGaussian());
        estimate = exact.getVariance();
      }
      else
      {
        LocalRFRunner rf = new LocalRFRunner(options);
        RealVariable monitored = spec.variables.get(0);
        rf.init(spec);
        rf.addMomentRayProcessor();
        rf.run();
        estimate = rf.momentRayProcessor.getVarianceEstimate(monitored);
      }
      mcErrorStatistics.addValue((estimate - 1.0)*(estimate - 1.0));
    }
    return Math.sqrt(mcErrorStatistics.getMean());
  }

  @Override
  public void run()
  {
    List<Double> xs = new ArrayList<>();
    List<Double> ys = new ArrayList<>();
    
    SimpleRegression reg = new SimpleRegression();
    for (int N = minN; N <= maxN; N *= 2)
    {
      double currentError = mcError(N);
      double logError = Math.log10(currentError);
      double logN = Math.log10(N);
      xs.add(logN);
      ys.add(logError);
      reg.addData(logN, logError);
    }
    output.printWrite("fit", "slope", reg.getSlope(), "intercept", reg.getIntercept(), "rSquare", reg.getRSquare());
    PlotLine.from(xs, ys).toPDF(Results.getFileInResultFolder("mcErrorScalings.pdf"));
    output.close();
  }

}
