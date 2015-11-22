package rejfree.processors;

import rejfree.local.LocalRFSampler;
import rejfree.local.TrajectoryRay;
import blang.variables.RealVariable;
import briefj.collections.Counter;



public class MomentRayProcessor implements RayProcessor
{
  private Counter<RealVariable> 
    sum   = new Counter<RealVariable>(),
    sumSq = new Counter<RealVariable>();
  private double currentTime = 0.0;

  public double getMeanEstimate(RealVariable variable) 
  { 
    return sum.getCount(variable) / currentTime; 
  }
  
  public double getVarianceEstimate(RealVariable variable)
  {
    final double muBar = getMeanEstimate(variable);
    return sumSq.getCount(variable) / currentTime - (muBar * muBar);
  }
  
  @Override
  public void init(LocalRFSampler sampler) {}
  
  @Override
  public void processRay(RealVariable var, TrajectoryRay ray, double time,
      LocalRFSampler sampler)
  {
    sum.incrementCount(var, 
        indefIntegralForMean(ray.position_t, ray.velocity_t, time - ray.t));
    
    sumSq.incrementCount(var, 
        indefIntegralForVar(ray.position_t, ray.velocity_t, time - ray.t));
    currentTime = time;
  }
  
  public static double indefIntegralForMean(double x0, double v, double t)
  {
    return x0 * t + v * t*t / 2.0;
  }
  
  public static double indefIntegralForVar(double x0, double v, double t)
  {
    return x0*x0 * t + x0 * v * t*t + v*v * t*t*t / 3.0;
  }
}