package rejfree.models.normal;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.SparseDoubleCholeskyDecomposition;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import rejfree.local.CollisionFactor;
import blang.annotations.DefineFactor;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.Factor;
import blang.factors.FactorList;
import blang.variables.RealVariable;
import briefj.BriefLists;



public class BrownianBridge
{
  public final List<Double> ts;
  public final List<RealVariable> variables;
  public final RealVariable globalVariance;
  
  public BrownianBridge(List<Double> ts, List<RealVariable> variables,
      RealVariable globalVariance)
  {
    this.ts = ts;
    this.variables = variables;
    this.globalVariance = globalVariance;
    globalVariance.setValue(1.0);
  }
  
  public static BrownianBridge regularlySpaced(int nVariables)
  {
    List<RealVariable> vars = new ArrayList<>();
    for (int i = 0; i < nVariables; i++)
      vars.add(new RealVariable(0.0));
    
    double delta = 1.0 / (1.0 + nVariables);
    List<Double> ts = new ArrayList<>();
    double cur = 0.0;
    for (int i = 0; i < nVariables + 1; i++)
    {
      cur += delta;
      ts.add(cur);
    }
    
    return new BrownianBridge(ts, vars, new RealVariable(1.0));
  }
  
  public double logAbsDetPrecision()
  {
    // compute Cholesky
    SparseDoubleCholeskyDecomposition chol = 
        new SparseDoubleCholeskyDecomposition(
            precisionMatrix().getColumnCompressed(false), 
            0); // use natural order
    DoubleMatrix2D L = chol.getL();
    
    double sum = 0.0;
    for (int i = 0; i < variables.size(); i++)
      sum += Math.log(Math.abs(L.get(i, i)));
    
    // *2 because det(precision) = det(L) * det(L^T)
    return 2*sum;
  }
  
  public SparseDoubleMatrix2D precisionMatrix()
  {
    final int nVariables = variables.size();
    SparseDoubleMatrix2D result = new SparseDoubleMatrix2D(nVariables, nVariables, nVariables * 3, 0.2, 0.5);
    
    // precision term for the first end of bridge
    result.set(0, 0, precision(ts.get(0)));
    
    // .. and the other 
    final int lastIndex = variables.size() - 1;
    result.set(lastIndex, lastIndex, precision(lastDelta()));
    
    // iterate over pairs of variable with index i and i+1
    for (int i = 0; i < variables.size() - 1; i++)
    {
      final double delta = ts.get(i + 1) - ts.get(i);
      double precision = precision(delta);
      increment(result, i + 0, i + 0, + precision); increment(result, i + 0, i + 1, - precision); 
      increment(result, i + 1, i + 0, - precision); increment(result, i + 1, i + 1, + precision); 
    }
    
    return result;
  }
  
  private static void increment(DoubleMatrix2D m, int row, int col, double increment)
  {
    m.set(row, col, increment + m.get(row, col));
  }
  
  public Factor fullFactor()
  {
    return new FullFactor(globalVariance, variables);
  }
  
  private class FullFactor implements Factor
  {
    @FactorArgument
    public final RealVariable globalVariable;
    
    @FactorComponent
    public final FactorList<RealVariable> normalMarginals;
    
    private FullFactor(RealVariable globalVariable,
        List<RealVariable> normalMarginals)
    {
      super();
      this.globalVariable = globalVariable;
      this.normalMarginals = FactorList.ofArguments(normalMarginals, true);
    }



    @Override
    public double logDensity()
    {
      // normalization
      
      double sum = 0.0;
      
      try 
      {
        sum = 0.5 * logAbsDetPrecision() - ((double) variables.size())/2.0 * Math.log(2.0 * Math.PI);
      }
      catch (IllegalArgumentException iae)
      {
        // if not pos def!
        return Double.NEGATIVE_INFINITY;
      }
      
      sum += quadraticForm(variables.get(0).getValue(), ts.get(0));
      sum += quadraticForm(BriefLists.last(variables).getValue(), lastDelta());
      
      for (int i = 0; i < variables.size() - 1; i++)
        sum += quadraticForm(
            variables.get(i).getValue() - variables.get(i + 1).getValue(), 
            ts.get(i + 1) - ts.get(i));
      
      return sum;
    }
  }
  
  private double quadraticForm(double x, double delta)
  {
    return -0.5 * x * x * precision(delta);
  }

  /**
   * Limitation: will become invalid if globalParam is updated.
   * 
   * For sampling globalParam, use fullFactor() instead.
   */
  public ArrayList<CollisionFactor> localFactors()
  {
    /*
     * Note 1: ideally may want to:
     * 
     *  1. construct a new NormalFactor where the precision is held in 
     *     a RealVariable instead AND
     *  2. modify LocalRFSampler to allow resuming iteration 
     *     with a recompute of collision times but no initial refreshment
     *  
     * But if it turns out that the rate at which we alternate between 
     * jump dynamics and collision dynamics is similar to the refresh rate,
     * might as well to a full re-initialization of collision time and 
     * velocity during the switch.
     */
    
    ArrayList<CollisionFactor> result = new ArrayList<>();
    
    // link to first end of the bridge
    result.add(buildUnary(ts.get(0) - 0.0, variables.get(0)));
    
    // link to second end of the bridge
    result.add(buildUnary(lastDelta(), BriefLists.last(variables)));
    
    // iterate over pairs of variable with index i and i+1
    for (int i = 0; i < variables.size() - 1; i++)
    {
      List<RealVariable> pairOfVariables = variables.subList(i, i + 2);
      
      final double delta = ts.get(i + 1) - ts.get(i);
      if (delta <= 0.0)
        throw new RuntimeException("Variables should be sorted and "
            + "have strictly increasing time indices.");
      final double precision = precision(delta);
      
      /*
       *  Note 2: not a true precisionMatrix (determinant = 0), but intermediate
       *        factors in L-BPS not required to be probability distributions 
       */
      DoubleMatrix pseudoPrecisionMatrix = new DoubleMatrix(2, 2,
          new double[]{
            + precision, - precision,
            - precision, + precision
          });
      
      result.add(new NormalFactor(
          pseudoPrecisionMatrix, 
          pairOfVariables, 
          0.0  // see Note 2 above
          ));
    }

    return result;
  }
  
  private double lastDelta()
  {
    return ts.get(ts.size() - 1) - ts.get(ts.size() - 2);
  }

  private CollisionFactor buildUnary(double delta, RealVariable var)
  {
    return NormalFactor.newUnaryFactor(new DoubleMatrix(new double[]{precision(delta)}), var);
  }
  
  private double precision(double delta)
  {
    return 1.0 / (delta * globalVariance.getValue());
  }
  
  public double getMarginalVariance(int i)
  {
    double t = ts.get(i);
    return globalVariance.getValue() * t * (1.0 - t);
  }

  public Object localFactorModelSpec()
  {
    return new Object() 
    {
      @DefineFactor
      public final List<CollisionFactor> factors = localFactors();
    };
  }
}
