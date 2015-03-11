package rejfree;

import org.jblas.Decompose;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import bayonet.opt.DifferentiableFunction;



public class NormalEnergy implements DifferentiableFunction
{
  private final DoubleMatrix precisionMatrix;
  
  /**
   * log((2pi)^{-k/2} + log(|sigma|^{-1/2}) 
   */
  private final double logConstant; 
  
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
    double detPrecAbs = Math.abs(Decompose.lu(precisionMatrix).u.diag().prod());
    this.logConstant = - (((double)precisionMatrix.getRows()) / 2.0) * Math.log(2.0 * Math.PI) + 0.5 * Math.log(detPrecAbs);
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
    return 0.5 * point.transpose().mmul(precisionMatrix).mmul(point).get(0) 
        - logConstant; // Note: this would be added to the log density, but is subtracted here because we need the energy
  }

  @Override
  public double[] derivativeAt(double[] x)
  {
    DoubleMatrix point = new DoubleMatrix(x);
    return precisionMatrix.mmul(point).data;
  }
  
}