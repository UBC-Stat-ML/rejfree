package rejfree;

import java.io.File;
import java.util.List;

import org.jblas.DoubleMatrix;

import bayonet.rplot.RJavaBridge;
import bayonet.rplot.RUtils;



public class PlotTrajectory extends RJavaBridge
{
  private String message = null;
  
  public double [] x0, y0, x1, y1;
  public double [] xBounds, yBounds;
  
  public static PlotTrajectory fromFirstTwoDimensions(List<DoubleMatrix> points)
  {
    return new PlotTrajectory(points, 0, 1);
  }
  
  public PlotTrajectory(List<DoubleMatrix> points, int xIndex, int yIndex)
  {
    double minX = Double.POSITIVE_INFINITY, maxX = Double.NEGATIVE_INFINITY, minY = Double.POSITIVE_INFINITY, maxY = Double.NEGATIVE_INFINITY;
    final int nPoints = points.size();
    x0 = new double[nPoints];
    y0 = new double[nPoints];
    x1 = new double[nPoints];
    y1 = new double[nPoints];
    for (int pointIndex = 0; pointIndex < points.size() - 1; pointIndex++)
    {
      DoubleMatrix 
        current = points.get(pointIndex),
        next    = points.get(pointIndex + 1);
      x0[pointIndex] = current.get(xIndex);  
      y0[pointIndex] = current.get(yIndex);
      x1[pointIndex] = next.get(xIndex);
      y1[pointIndex] = next.get(yIndex);
      if (x0[pointIndex] < minX)
        minX = x0[pointIndex];
      if (x0[pointIndex] > maxX)
        maxX = x0[pointIndex];
      if (y0[pointIndex] < minY)
        minY = y0[pointIndex];
      if (y0[pointIndex] > maxY)
        maxY = y0[pointIndex];
    }
    
    xBounds = new double[]{minX, maxX};
    yBounds = new double[]{minY, maxY};
  }

  /**
   * Save to pdf 
   * 
   * @param output
   */
  public void toPDF(File output)
  {
    this.output = output;
    message  = RUtils.callRBridge(this);
  }

  /**
   * @return The message returned by last call to R
   */
  public String getMessage()
  {
    return message;
  }
 
  @Override public String rTemplateResourceURL() { return "/rejfree/PlotTrajectory.txt"; }

  private File output;
  public String getOutput()
  {
    return RUtils.escapeQuote(output.getAbsolutePath());
  }
}
