package rejfree.local;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import rejfree.models.normal.IsotropicNormal;
import bayonet.coda.EffectiveSize;
import bayonet.rplot.PlotHistogram;
import bayonet.rplot.PlotLine;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;



public class CompareESSLocalGlobal implements Runnable
{
  @Option
  public boolean useSparse = false;
  
  @Option
  public double samplingInterval = 1;
  
  @Option
  public int minDim = 8;
  
  @Option
  public int maxDim = 1024;
  
  @Option
  public int nRepeatsForSlopeEstimation = 2;
  
  @Option
  public int nRepeatsForSlopeEstimatorVariability = 100;
  
  @Option
  public int nBurnIn = 1;

  @OptionSet(name = "rf")
  public LocalRFRunnerOptions options = new LocalRFRunnerOptions();
  
  @Override
  public void run()
  {
    List<Double> estimatedSlopes = new ArrayList<>();
    OutputManager out = Results.getGlobalOutputManager();
    for (int j = 0; j < nRepeatsForSlopeEstimatorVariability; j++)
    {
      List<Double> xs = new ArrayList<>();
      List<Double> ys = new ArrayList<>();
      
      SimpleRegression reg = new SimpleRegression();
      for (int dim = minDim; dim <= maxDim; dim *= 2)
      {
        SummaryStatistics stat = new SummaryStatistics();
        ModelSpec spec = new ModelSpec(dim, useSparse);
        for (int i = 0; i < nRepeatsForSlopeEstimation; i++)
        {
          LocalRFRunner rf = new LocalRFRunner(options);
          RealVariable monitored = spec.variables.get(0);
          rf.init(spec);
          rf.addSaveRaysProcessor(Collections.singleton(monitored));
          rf.run();
          List<Double> convertToSample = rf.saveRaysProcessor.convertToSample(monitored, samplingInterval);
          double ess = EffectiveSize.effectiveSize(convertToSample);
          long timems = rf.watch.elapsed(TimeUnit.MILLISECONDS);
          double essPerSec = (1000.0*ess/timems);
          out.printWrite("essPerSec", "dim", dim, "innerRep", i, "outerRep", j, "essPerSec", essPerSec, "ess", ess, "timeSec", (timems/1000.0));
          if (i >= nBurnIn)
            stat.addValue(essPerSec);
        }
        double meanEss = stat.getMean();
        reg.addData(Math.log10(dim), Math.log10(meanEss));
        xs.add(Math.log10(dim));
        ys.add(Math.log10(meanEss));
      }
      out.printWrite("fit", "outerRep", j, "slope", reg.getSlope(), "intercept", reg.getIntercept(), "rSquare", reg.getRSquare());
      PlotLine.from(xs, ys).toPDF(Results.getFileInResultFolder("essPerSecByDim-" + j + ".pdf"));
      estimatedSlopes.add(reg.getSlope());
      if (estimatedSlopes.size() > 1)
        PlotHistogram.from(estimatedSlopes).toPDF(Results.getFileInResultFolder("estimatedSlopes-estimatorSamplingDistribution.pdf"));
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
