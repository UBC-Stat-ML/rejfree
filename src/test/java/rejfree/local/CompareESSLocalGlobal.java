package rejfree.local;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import rejfree.models.normal.IsotropicNormal;
import bayonet.coda.EffectiveSize;
import bayonet.rplot.PlotLine;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class CompareESSLocalGlobal implements Runnable
{
  @Option
  public boolean useSparse = false;
  
  @Option
  public double samplingInterval = 0.1;

  
  @Override
  public void run()
  {
    LocalRFRunnerOptions options = new LocalRFRunnerOptions();
    options.maxRunningTimeMilli = 5*1000;
    options.maxSteps = Integer.MAX_VALUE;
    options.maxTrajectoryLength = Double.POSITIVE_INFINITY;
    
    List<Double> xs = new ArrayList<>();
    List<Double> ys = new ArrayList<>();
    
    SimpleRegression reg = new SimpleRegression();
    for (int dim = 32; dim < 1000000; dim *= 2)
    {
      SummaryStatistics stat = new SummaryStatistics();
      ModelSpec spec = new ModelSpec(dim, useSparse);
      for (int i = 0; i < 100; i++)
      {
        LocalRFRunner rf = new LocalRFRunner(options);
        RealVariable monitored = spec.variables.get(0);
        rf.init(spec);
        rf.addSaveRaysProcessor(Collections.singleton(monitored));
        rf.run();
        List<Double> convertToSample = rf.saveRaysProcessor.convertToSample(monitored, samplingInterval);
        double ess = EffectiveSize.effectiveSize(convertToSample);
        long timems = rf.watch.elapsed(TimeUnit.MILLISECONDS);
        rf.output.printWrite("ess", "ess", ess);
        double essPerSec = (1000.0*ess/timems);
        rf.output.printWrite("essPerSec", "dim", dim, "essPerSec", essPerSec, "ess", ess, "timeSec", (timems/1000.0));
        if (i > 3)
          stat.addValue(essPerSec);
      }
      double meanEss = stat.getMean();
      reg.addData(Math.log10(dim), Math.log10(meanEss));
      xs.add(Math.log10(dim));
      ys.add(Math.log10(meanEss));
      System.out.println("slope=" + reg.getSlope());
      PlotLine.from(xs, ys).toPDF(Results.getFileInResultFolder("essPerSecByDim.pdf"));
    }
  }
  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new CompareESSLocalGlobal());
  }
  
  public class ModelSpec
  {
    public final List<RealVariable> variables = new ArrayList<RealVariable>();
    
    @DefineFactor
    public final List<CollisionFactor> factors = new ArrayList<>();
    
    public ModelSpec(int dim, boolean sparse)
    {
      for (int i = 0; i < dim; i++)
        variables.add(new RealVariable(0.0));
      if (sparse)
      {
        for (int i = 0; i < dim; i++)
          factors.add(new IsotropicNormal(Collections.singletonList(variables.get(i))));
      }
      else
        this.factors.add(new IsotropicNormal(variables));
    }
  }
}
