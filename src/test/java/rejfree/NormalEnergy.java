package rejfree;

import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import bayonet.opt.DifferentiableFunction;



public class NormalEnergy implements DifferentiableFunction
{
  private final DoubleMatrix precisionMatrix;
  
  public static NormalEnergy isotropic(int dim)
  {
    DoubleMatrix precision = DoubleMatrix.eye(dim);
    return new NormalEnergy(precision);
  }
  
  public static NormalEnergy withCovariance(DoubleMatrix covar)
  {
    return new NormalEnergy(Solve.pinv(covar));
  }
  
  private NormalEnergy(DoubleMatrix precisionMatrix)
  {
    this.precisionMatrix = precisionMatrix;
  }

  @Override
  public int dimension()
  {
    return precisionMatrix.getRows();
  }

  @Override
  public double valueAt(double[] x)
  {
    DoubleMatrix point = new DoubleMatrix(x);
    return 0.5 * point.transpose().mmul(precisionMatrix).mmul(point).get(0);
  }

  @Override
  public double[] derivativeAt(double[] x)
  {
    DoubleMatrix point = new DoubleMatrix(x);
    return precisionMatrix.mmul(point).data;
  }
  
}