package rejfree;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class HMCTest implements Runnable
{
  @Option
  public static double epsilon = 0.3;
  
  @Option
  public static int nIters = 100;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new HMCTest());
  }
  
  public static List<DoubleMatrix> doHMC(Function<Pair<DoubleMatrix,DoubleMatrix>,Pair<DoubleMatrix,DoubleMatrix>> integrator)
  {
    DoubleMatrix 
    q = new DoubleMatrix(new double[]{0,1}),
    p = new DoubleMatrix(new double[]{1,0});
    
    Pair<DoubleMatrix,DoubleMatrix> current = Pair.of(q, p);
    List<DoubleMatrix> trajectory = new ArrayList<>();
    trajectory.add(current.getLeft());
    
    for (int i = 0; i < nIters; i++)
    {
      current = integrator.apply(current);
      trajectory.add(current.getLeft());
    }
    
    return trajectory;
  }
  
  public static Pair<DoubleMatrix,DoubleMatrix> leapFrog(Pair<DoubleMatrix,DoubleMatrix> current)
  {
    DoubleMatrix 
      q = current.getLeft().dup(),
      p = current.getRight().dup();
    
    p.subi(q.mul(epsilon/2.0));
    q.addi(p.mul(epsilon));
    p.subi(q.mul(epsilon/2.0));
    
    return Pair.of(q, p);
  }
  
  public static Pair<DoubleMatrix,DoubleMatrix> euler(Pair<DoubleMatrix,DoubleMatrix> current)
  {
    DoubleMatrix 
      q = current.getLeft().dup(),
      p = current.getRight().dup();
    
    p.subi(q.mul(epsilon));
    q.addi(p.mul(epsilon));
    
    return Pair.of(q, p);
  }

  @Override
  public void run()
  {
    PlotTrajectory.fromFirstTwoDimensions(doHMC(HMCTest::euler)).toPDF(Results.getFileInResultFolder("euler.pdf"));
    PlotTrajectory.fromFirstTwoDimensions(doHMC(HMCTest::leapFrog)).toPDF(Results.getFileInResultFolder("leapFrog.pdf"));
  }

}
