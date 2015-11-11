package rejfree.mixture;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.analysis.MultivariateFunction;

import rejfree.models.normal.NormalFactor;
import bayonet.rplot.PlotContour;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
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
    public final NormalFactor priorForComponentMean0 = NormalFactor.newUnaryStandardNormal(meanForComponentMean0);
    
    @DefineFactor
    public final NormalFactor priorForComponentMean1 = NormalFactor.newUnaryStandardNormal(meanForComponentMean1);
    
    @DefineFactor(onObservations = true)
    public final List<MixtureFactor> likelihood;
    
    private MixtureSpec()
    {
      likelihood = new ArrayList<MixtureFactor>(nObservations);
      for (int i = 0; i < nObservations; i++)
      {
        MixtureFactor mixtureFactor = new MixtureFactor(meanForComponentMean0, meanForComponentMean1);
        likelihood.add(mixtureFactor);
        mixtureFactor.generate(random);
      }
    }
    
    
  }
}
