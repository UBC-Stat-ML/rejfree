package rejfree.local;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jblas.DoubleMatrix;

import com.google.common.base.Stopwatch;

import rejfree.RFSamplerOptions;
import rejfree.RFSamplerOptions.RefreshmentMethod;
import rejfree.global.GlobalRFSampler;
import rejfree.global.GlobalRFSampler.CollisionSolver;
import rejfree.models.normal.IsotropicNormal;
import rejfree.models.normal.NormalFactor;
import bayonet.coda.EffectiveSize;
import bayonet.opt.DifferentiableFunction;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.opt.Option;
import briefj.run.Mains;



public class CompareESSLocalGlobal implements Runnable
{
  
  public static class NormalCollisionSolver implements CollisionSolver
  {

    @Override
    public double collisionTime(DoubleMatrix x,
        DoubleMatrix v, DifferentiableFunction energy, double e)
    {
      final double xv = x.dot(v);
      final double vv = v.dot(v);
      return NormalFactor.normalCollisionTime(e, xv, vv);
    }
    
  }
  
  @Override
  public void run()
  {
//    final int dim = 10;
//    DoubleMatrix initialPosition = new DoubleMatrix(dim);
//    RFSamplerOptions options = new RFSamplerOptions();
//    options.refreshmentMethod = RefreshmentMethod.GLOBAL;
//    
//    CollisionSolver solver = new NormalCollisionSolver();
//    GlobalRFSampler sampler = new GlobalRFSampler(new DifferentiableFunction() {
//      
//      @Override
//      public double valueAt(double[] x)
//      {
//        double sum = 0.0;
//        for (int i = 0; i < x.length; i++)
//          sum += x[i] * x[i];
//        return 0.5 * sum;
//      }
//      
//      @Override
//      public int dimension()
//      {
//        return dim;
//      }
//      
//      @Override
//      public double[] derivativeAt(double[] x)
//      {
//        return x;
//      }
//    }, initialPosition , options , solver );
//    
//    for (int rep = 0; rep < 10; rep++)
//    {
//      Random rand = new Random(1);
//      Stopwatch watch = Stopwatch.createStarted();
//      sampler.iterate(rand , 10000);
//      long elapsed = watch.elapsed(TimeUnit.MILLISECONDS);
//      List<DoubleMatrix> samples = sampler.getSamples();
//      List<Double> firsts = new ArrayList<>();
//      for (int i = 0; i < samples.size(); i++)
//        firsts.add(samples.get(i).get(0));
//      double ess = EffectiveSize.effectiveSize(firsts);
//      System.out.println("elapsed(ms)=" + elapsed);
//      System.out.println("ess=" + ess);
//      System.out.println("ess/ms=" + (ess/elapsed));
//    }
    
    SimpleRegression reg = new SimpleRegression();
    for (int dim = 10; dim < 10000; dim *= 10)
    {
      SummaryStatistics stat = new SummaryStatistics();
      for (int i = 0; i < 10; i++)
      {
        LocalRFRunner rf = new LocalRFRunner();
        rf.options.maxRunningTimeMilli = 500;
        rf.options.maxSteps = Integer.MAX_VALUE;
        rf.options.maxTrajectoryLength = Double.POSITIVE_INFINITY;
        ModelSpec spec = new ModelSpec(dim, false);
        RealVariable monitored = spec.variables.get(0);
        rf.init(spec);
        rf.addSaveRaysProcessor(Collections.singleton(monitored));
        rf.run();
        List<Double> convertToSample = rf.saveRaysProcessor.convertToSample(monitored, 0.1);
        double ess = EffectiveSize.effectiveSize(convertToSample);
        long timems = rf.watch.elapsed(TimeUnit.MILLISECONDS);
        rf.output.printWrite("ess", "ess", ess);
        double essPerSec = (1000.0*ess/timems);
        rf.output.printWrite("essPerSec", "essPerSec", essPerSec);
        if (i > 3)
          stat.addValue(essPerSec);
      }
      double meanEss = stat.getMean();
      reg.addData(Math.log(dim), Math.log(meanEss));
      System.out.println("slope=" + reg.getSlope());
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
