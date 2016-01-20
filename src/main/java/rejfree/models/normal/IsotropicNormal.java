package rejfree.models.normal;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import rejfree.StaticUtils;
import rejfree.local.CollisionContext;
import rejfree.local.CollisionFactor;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.factors.GenerativeFactor;
import blang.variables.RealVariable;


/**
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class IsotropicNormal implements CollisionFactor, GenerativeFactor
{
  @FactorComponent
  public final FactorList<RealVariable> variables;
  
  private final int dim;
  
  public IsotropicNormal(List<RealVariable> variables)
  {
    this.dim = variables.size();
    this.variables = FactorList.ofArguments(variables, true);
  }

  @Override
  public double logDensity()
  {
    DoubleMatrix point = getPosition();
    return - 0.5 * dotProd(point, point); 
  }
  
  private double dotProd(final DoubleMatrix x1, final DoubleMatrix x2)
  {
    return x1.dot(x2);
  }
  
  public int dim()
  {
    return dim;
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
    
    return Pair.of(NormalFactor.normalCollisionTime(e, xv, vv),true); 
  }
  


  @Override
  public DoubleMatrix gradient()
  {
    final DoubleMatrix x = getPosition();
    return x.mul(-1.0);
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
    for (int i = 0; i < dim; i++)
      variables.list.get(i).setValue(random.nextGaussian());
  }
}
