package rejfree;

import blang.variables.RealVariable;



public interface RayProcessor
{
  public void init(LocalRFSampler sampler);
  public void processRay(RealVariable var, TrajectoryRay ray, double timeTheRayEnds, LocalRFSampler sampler);
}