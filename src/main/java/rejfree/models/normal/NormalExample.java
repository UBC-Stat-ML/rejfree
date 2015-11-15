package rejfree.models.normal;

import java.io.File;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;

import rejfree.PlotTrajectory;
import rejfree.global.GlobalRFSampler;
import rejfree.global.GlobalRFSampler.RFSamplerOptions;

import com.google.common.collect.Lists;

import bayonet.coda.EffectiveSize;
import bayonet.rplot.PlotHistogram;
import blang.MCMCFactory.MCMCOptions;
import blang.processing.ProcessorContext;
import blang.variables.RealVariable;
import blang.variables.RealVariableProcessor;
import briefj.BriefStrings;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;

public class NormalExample implements Runnable
{
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new NormalExample());
  }

  @Option
  public int nIterations = 10000;
  
  @Option
  public int maxTrajectoryToPlot = 500;
  
  @OptionSet(name = "rf")
  public RFSamplerOptions samplerOptions = new RFSamplerOptions();
  
  @Option
  public double var1 = 1.0;
  
  @Option
  public double var2 = 1.0;
  
  @Option
  public double cov = 0.5;

  @Override
  public void run()
  {
    DoubleMatrix covar = new DoubleMatrix(new double[][]{{var1,cov},{cov,var2}});
    NormalEnergy energy = NormalEnergy.
      withCovariance(covar); 
      //isotropic(2);
    
    GlobalRFSampler sampler = GlobalRFSampler.initializeRFWithLBFGS(energy, samplerOptions);
    Random rand = new Random(134);
    System.out.println("Sampling " + nIterations + " events");
    sampler.iterate(rand, nIterations);
    
    // print some statistics
    if (maxTrajectoryToPlot > 2)
    {
      File output = Results.getFileInResultFolder("trajectory.pdf");
      System.out.println("Printing trajectory of " + maxTrajectoryToPlot + " first events in " + output);
      PlotTrajectory pt = new PlotTrajectory(sampler.getTrajectory().subList(0, Math.min(maxTrajectoryToPlot, sampler.getTrajectory().size())), 0, 1);
      pt.toPDF(output);
    }
    
    // monitoring behavior of event sampler
    System.out.println("Mean collision to refreshment ratio: " + sampler.getCollisionToRefreshmentRatio().getMean());
    System.out.println("Mean # collected events per event: " + sampler.getCollectedPerEvent().getMean());
    System.out.println("Number of samples collected: " + sampler.getSamples().size());
    
    // estimates
    System.out.println();
    DoubleMatrix meanEstimate = new DoubleMatrix(sampler.dimensionality());
    for (DoubleMatrix sample : sampler.getSamples())
      meanEstimate.addi(sample);
    meanEstimate.muli(1.0/sampler.getSamples().size());
    System.out.println("MC estimate of mean:\n" + BriefStrings.indent(meanEstimate.toString()));
    DoubleMatrix covarEstimate = DoubleMatrix.zeros(sampler.dimensionality(), sampler.dimensionality());
    for (DoubleMatrix sample : sampler.getSamples())
    {
      DoubleMatrix centered = sample.sub(meanEstimate);
      covarEstimate.addi(centered.mmul(centered.transpose()));
    }
    covarEstimate.muli(1.0/sampler.getSamples().size());
    System.out.println("MC estimate of covar:\n" + BriefStrings.indent(covarEstimate.toString()));
    
    // print info on marginals
    System.out.println();
    File marginalHistDir = Results.getFolderInResultFolder("marginal-histograms");
    System.out.println("Printing marginal histograms in " + marginalHistDir);
    System.out.println("mean from full traj: " + sampler.getMean());
    System.out.println("var from full traj: " + sampler.getVariance());
    for (int d = 0; d < 2; d++)
    {
      List<Double> coordinates = Lists.newArrayList();
      SummaryStatistics marginalStat = new SummaryStatistics();
      RealVariable dummy = RealVariable.real();
      RealVariableProcessor processor = new RealVariableProcessor("marginal-" + d, dummy);
      int iter = 0;
      MCMCOptions options = new MCMCOptions();
      options.nMCMCSweeps = sampler.getSamples().size();
      for (DoubleMatrix sample : sampler.getSamples())
      {
        final double current = sample.get(d);
        marginalStat.addValue(current);
        coordinates.add(current);
        dummy.setValue(current);
        
        ProcessorContext context = new ProcessorContext(iter++, null, options);
        processor.process(context);
      }
      PlotHistogram.from(coordinates).toPDF(new File(marginalHistDir, "hist-" + d + ".pdf"));
      System.out.println("Marginal " + d);
      System.out.println(BriefStrings.indent(marginalStat.toString()));
      System.out.println("ESS: " + EffectiveSize.effectiveSize(coordinates));
      System.out.println();
    }
  }
}
