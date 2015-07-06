package rejfree.glm;

import org.jblas.DoubleMatrix;

import bayonet.math.JBlasUtils;



public class TestMatrix
{

  public static void main(String[] args)
  {
    final double a = 0.5, b = 0.4, c = 0.2;
    DoubleMatrix m = new DoubleMatrix(new double[][]{{1,a,0,0},{a,1,b,0},{0,b,1,c},{0,0,c,1}});
    
    System.out.println(JBlasUtils.inversePositiveMatrix(m));
    
    
  }

}
