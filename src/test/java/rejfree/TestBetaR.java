package rejfree;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;

import rejfree.GlobalRFSampler.RFSamplerOptions;




public class TestBetaR
{
  private static class Event
  {
    double t, beta, R;
  }
  
  public static List<Event> convertToEvents(List<DoubleMatrix> trajectory)
  {
    double timeConsumed = 0.0;
    List<Event> events = new ArrayList<Event>();
    for (int i = 1; i < trajectory.size() - 1; i++)
    {
      DoubleMatrix 
          prev = trajectory.get(i-1),
          cur  = trajectory.get(i),
          next = trajectory.get(i+1);
     
      DoubleMatrix 
        diffAfter  = next.sub(cur),
        diffBefore = cur.sub(prev);
//      DoubleMatrix veloBefore = diffBefore.div(diffBefore.norm2());
      DoubleMatrix veloAfter  = diffAfter.div(diffAfter.norm2());
      timeConsumed += diffBefore.norm2();
      double beta = Math.acos(cur.dot(veloAfter) / cur.norm2());
      
      Event e = new Event();
      e.beta = beta;
      e.t = timeConsumed;
      e.R = cur.norm2();
      
      events.add(e);
    }
    return events;
  }
  
  public static double computeRPrime(double R0, double beta0, double t)
  {
    return Math.sqrt( R0 * R0 + 2 * t * Math.cos(beta0) * R0 + t * t);
  }
  
  public static double computeBetaPrime(double R0, double beta0, double t)
  {
    return Math.PI - Math.acos(  (Math.cos(beta0) * R0 + t)  /  computeRPrime(R0, beta0, t) );
  }
  
  public static void main(String [] args)
  {
    final int dim = 2;
    final int nCollisions = 1000;
    
    RFSamplerOptions options = new RFSamplerOptions();
    Random rand = new Random(1);
    NormalEnergy e = NormalEnergy.isotropic(dim);
    DoubleMatrix initPos = DoubleMatrix.randn(dim);
    options.refreshRate = 0.0;
    options.collectRate = 0.0;
    GlobalRFSampler sampler = new GlobalRFSampler(e, initPos.dup() , options);
      
    sampler.iterate(rand , nCollisions);
    List<Event> trajectory = convertToEvents(sampler.getTrajectory());
    
    for (int i = 1; i < trajectory.size() - 1; i++)
    {
      System.out.println("Iteration " + i);
      Event 
        current = trajectory.get(i),
        next    = trajectory.get(i+1);
      
      double analyticRPrime = computeRPrime(current.R, current.beta, next.t - current.t);
      double observedRPrime = next.R;
      
      double analyticBetaPrime = computeBetaPrime(current.R, current.beta, next.t - current.t);
      double observedBetaPrime = next.beta;
      
      if (analyticBetaPrime < Math.PI /2.0 || observedBetaPrime < Math.PI / 2.0)
        System.err.println("Warning");
      
      System.out.println("\tR:");
      System.out.println("\t\tAnalytic: " + analyticRPrime);
      System.out.println("\t\tObserved: " + observedRPrime);
      
      System.out.println("\tbeta:");
      System.out.println("\t\tAnalytic: " + analyticBetaPrime);
      System.out.println("\t\tObserved: " + observedBetaPrime);
    }
  }
}
