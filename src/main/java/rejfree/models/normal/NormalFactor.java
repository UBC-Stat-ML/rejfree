package rejfree.models.normal;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;

import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
import bayonet.math.JBlasUtils;
import bayonet.math.NumericalUtils;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.factors.GenerativeFactor;
import blang.variables.RealVariable;


/**
 * NOTE: this is used in the local sampler, so we do NOT assume velocity for 
 *   the variables of interest to be of unit norm
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class NormalFactor implements CollisionFactor, GenerativeFactor
{
  @FactorComponent
  public final FactorList<RealVariable> variables;
  
  // Note: this implementation would not currently support making precision random
  private final DoubleMatrix precision;
  
  private final boolean isBin;
  private final double p0, p1, d;
  
  /**
   * log((2pi)^{-k/2} + log(|sigma|^{-1/2}) 
   */
  private final double logConstant; 
  
  public static NormalFactor newUnaryFactor(DoubleMatrix precision, RealVariable variable0)
  {
    List<RealVariable> variables = new ArrayList<RealVariable>();
    variables.add(variable0);
    return new NormalFactor(precision, variables);
  }
  
  public static NormalFactor newBinaryFactor(DoubleMatrix precision, RealVariable variable0, RealVariable variable1)
  {
    List<RealVariable> variables = new ArrayList<RealVariable>();
    variables.add(variable0);
    variables.add(variable1);
    return new NormalFactor(precision, variables);
  }
  
  public static NormalFactor withCovariance(DoubleMatrix covar, List<RealVariable> variables)
  {
    return new NormalFactor(JBlasUtils.inversePositiveMatrix(covar), variables);
  }
  
  public NormalFactor(DoubleMatrix precision, List<RealVariable> variables)
  {
    this.variables = FactorList.ofArguments(variables, true);
    this.precision = precision.dup();
    double detPrecAbs = Math.abs(Decompose.lu(precision).u.diag().prod());
    this.logConstant = - (((double)precision.getRows()) / 2.0) * Math.log(2.0 * Math.PI) + 0.5 * Math.log(detPrecAbs);
    isBin = (variables.size() == 2);
    p0 = isBin ? precision.get(0,0) : Double.NaN;
    p1 = isBin ? precision.get(1,1) : Double.NaN;
    d  = isBin ? precision.get(0,1) : Double.NaN;
    if (isBin)
      NumericalUtils.checkIsClose(precision.get(0,1), precision.get(1,0));
  }

  @Override
  public double logDensity()
  {
    DoubleMatrix point = getPosition();
    return - 0.5 * dotProd(point, point) + logConstant; 
  }
  
  private double dotProd(final DoubleMatrix x1, final DoubleMatrix x2)
  {
    if (isBin)
    { // optimization of computational inner-loop bottleneck
      final double [] 
        array0 = x1.data,
        array1 = x2.data;
      return 
          array0[0] * (array1[0] * p0 + array1[1] * d) + 
          array0[1] * (array1[0] * d  + array1[1] * p1);
//  // Test of correctness:
//      double v1 = array0[0] * (array1[0] * p0 + array1[1] * d) + 
//             array0[1] * (array1[0] * d  + array1[1] * p1);
//      double v2 = x1.transpose().mmuli(precision).dot(x2);
//      NumericalUtils.checkIsClose(v1, v2);
//      if (v1 != v2)
//        System.err.println("slight difference: " + v1 + " vs " + v2 + " delta = " + (v1 - v2));
//      return v2;
    }
    else
      return x1.transpose().mmuli(precision).dot(x2);
  }
  
  public int dim()
  {
    return precision.rows;
  }

  private DoubleMatrix getPosition()
  {
    DoubleMatrix result = new DoubleMatrix(dim());
    for (int i = 0; i < dim(); i++)
      result.put(i, variables.list.get(i).getValue());
    return result;
  }
  
  public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(
      CollisionContext context)
  {
    final DoubleMatrix x = getPosition();
    final DoubleMatrix v = context.velocity;
    
    final double xv = dotProd(x, v);
    final double vv = dotProd(v, v);
    final double e = StaticUtils.generateUnitRateExponential(context.random);
    
    return Pair.of(normalCollisionTime(e, xv, vv),true); 
  }
  
  public static double normalCollisionTime(double exponential, double xv, double vv)
  {
    final double s1 = xv < 0 ? - xv / vv : 0.0;
    final double C = - exponential - s1 * (xv + vv * s1 / 2.0);
    final double result = (- xv + Math.sqrt(xv * xv - 2.0 * vv * C)) / vv;
    return result;
  }

  @Override
  public DoubleMatrix gradient()
  {
    final DoubleMatrix x = getPosition();
    return precision.mmul(x).mul(-1.0);
  }

  @Override
  public RealVariable getVariable(int gradientCoordinate)
  {
    return variables.list.get(gradientCoordinate);
  }

  @Override
  public int nVariables()
  {
    return dim();
  }

  @Override
  public void generate(Random random)
  {
    DoubleMatrix covarMatrix = JBlasUtils.inversePositiveMatrix(precision);  
    MultivariateNormalDistribution normal = new MultivariateNormalDistribution(new DoubleMatrix(dim()).data, JBlasUtils.asDoubleArray(covarMatrix));
    normal.reseedRandomGenerator(random.nextLong());
    double[] sample = normal.sample();
    for (int i = 0; i < sample.length; i++)
      variables.list.get(i).setValue(sample[i]);
  }
}
