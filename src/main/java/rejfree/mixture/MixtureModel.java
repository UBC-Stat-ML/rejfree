package rejfree.mixture;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.jblas.DoubleMatrix;

import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
import rejfree.models.normal.NormalFactor;
import bayonet.rplot.PlotContour;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.annotations.FactorArgument;
import blang.variables.RealVariable;
import briefj.opt.Option;


/**
 * Warning: this was not implemented with the goal of efficiency or
 * generality in mind, but rather to quickly create intuitive plots
 * explaining the execution of the algorithm.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class MixtureModel
{
  @Option
  public Random random = new Random(1);

  @Option
  public int nObservations = 20;
  
  @Option
  public double trueMean0 = -1;
  
  @Option
  public double trueMean1 = +1;
  
  @Option
  public boolean useLocal = false;
  
  public MixtureSpec getModelSpec()
  {
    return new MixtureSpec();
  }
  
  public PlotContour densityPlot()
  {
    MixtureSpec modelSpec = getModelSpec();
    ProbabilityModel m = new ProbabilityModel(modelSpec);
    
    PlotContour plot = new PlotContour(new MultivariateFunction() {
      
      @Override
      public double value(double[] point)
      {
        modelSpec.meanForComponentMean0.setValue(point[0]);
        modelSpec.meanForComponentMean1.setValue(point[1]);
        return m.logDensity();
      }
    });
    final double 
      min = Math.min(trueMean0, trueMean1) - 1,
      max = Math.max(trueMean0, trueMean1) + 1;
    
    plot.min_x = min;
    plot.min_y = min;
    plot.max_x = max;
    plot.max_y = max;
    
    return plot;
  }
  
  public class MixtureSpec
  {
    public final RealVariable meanForComponentMean0 = RealVariable.real(trueMean0);
    public final RealVariable meanForComponentMean1 = RealVariable.real(trueMean1);
    
    @DefineFactor
    public NormalFactor priorForComponentMean0;
    
    @DefineFactor
    public NormalFactor priorForComponentMean1;
    
    @DefineFactor(onObservations = true)
    public List<MixtureFactor> likelihood;
    
    @DefineFactor
    public CombinedFactor combinedFactor;
    
    private MixtureSpec()
    {
      priorForComponentMean0 = NormalFactor.newUnaryStandardNormal(meanForComponentMean0);
      priorForComponentMean1 = NormalFactor.newUnaryStandardNormal(meanForComponentMean1);
      likelihood = new ArrayList<MixtureFactor>(nObservations);
      for (int i = 0; i < nObservations; i++)
      {
        MixtureFactor mixtureFactor = new MixtureFactor(meanForComponentMean0, meanForComponentMean1);
        likelihood.add(mixtureFactor);
        mixtureFactor.generate(random);
      }
      if (useLocal)
        combinedFactor = null;
      else
      {
        combinedFactor = new CombinedFactor(priorForComponentMean0, priorForComponentMean1, likelihood);
        priorForComponentMean0 = null;
        priorForComponentMean1 = null;
        likelihood = null;
      }
    }
  }
  
  public static class CombinedFactor implements CollisionFactor
  {
    private final List<MixtureFactor> likelihood;
    private final NormalFactor priorForComponentMean1;
    private final NormalFactor priorForComponentMean0;
    
    @FactorArgument(makeStochastic = true)
    public final RealVariable mean0;
    
    @FactorArgument(makeStochastic = true)
    public final RealVariable mean1;

    public CombinedFactor(NormalFactor priorForComponentMean0,
        NormalFactor priorForComponentMean1, List<MixtureFactor> likelihood)
    {
      this.priorForComponentMean0 = priorForComponentMean0;
      this.priorForComponentMean1 = priorForComponentMean1;
      this.likelihood = likelihood;
      this.mean0 = priorForComponentMean0.getVariable(0);
      this.mean1 = priorForComponentMean1.getVariable(0); // 0 correct in last clause
    }

    @Override
    public double logDensity()
    {
      double sum = 0.0;
      sum += priorForComponentMean0.logDensity();
      sum += priorForComponentMean1.logDensity();
      for (MixtureFactor f : likelihood)
        sum += f.logDensity();
      return sum;
    }

    @Override
    public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
        CollisionContext context)
    {
      double min = Double.POSITIVE_INFINITY;
      Pair<Double,Boolean> argmin = null;
      
      {
        Pair<Double, Boolean> lowerBoundForCollisionDeltaTime0 = priorForComponentMean0.getLowerBoundForCollisionDeltaTime(context);
        if (lowerBoundForCollisionDeltaTime0.getLeft() < min)
        {
          argmin = lowerBoundForCollisionDeltaTime0;
          min = lowerBoundForCollisionDeltaTime0.getLeft();
        }
      }
      
      {
        Pair<Double, Boolean> lowerBoundForCollisionDeltaTime1 = priorForComponentMean1.getLowerBoundForCollisionDeltaTime(context);
        if (lowerBoundForCollisionDeltaTime1.getLeft() < min)
        {
          argmin = lowerBoundForCollisionDeltaTime1;
          min = lowerBoundForCollisionDeltaTime1.getLeft();
        }
      }
      
      for (MixtureFactor f : likelihood)
      {
        Pair<Double, Boolean> item = f.getLowerBoundForCollisionDeltaTime(context);
        if (item.getLeft() < min)
        {
          argmin = item;
          min = item.getLeft();
        }
      }
      
      return argmin;
    }

    @Override
    public DoubleMatrix gradient()
    {
      DoubleMatrix result = new DoubleMatrix(2);
      
      result.put(0, priorForComponentMean0.gradient().get(0));
      result.put(1, priorForComponentMean1.gradient().get(0)); // 0 is correct in last clause
      
      for (MixtureFactor f : likelihood)
        result.addi(f.gradient());
      
      return result;
    }

    @Override
    public RealVariable getVariable(int gradientCoordinate)
    {
      return gradientCoordinate == 0 ? priorForComponentMean0.getVariable(0) : priorForComponentMean1.getVariable(0); // 0 is correct in second clause
    }

    @Override
    public int nVariables()
    {
      return 2;
    }
    
  }
}
