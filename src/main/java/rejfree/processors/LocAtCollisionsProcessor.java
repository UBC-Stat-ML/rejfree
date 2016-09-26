package rejfree.processors;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;

import org.jblas.DoubleMatrix;

import rejfree.local.LocalRFSampler;
import rejfree.local.TrajectoryRay;
import blang.variables.RealVariable;



public class LocAtCollisionsProcessor implements RayProcessor
{
  public final LinkedHashSet<RealVariable> variables;
  public final List<DoubleMatrix> locationsAtCollisions = new ArrayList<>();
  private double prevTime = -1;
  private LocalRFSampler sampler = null;

  public LocAtCollisionsProcessor(RealVariable ... variables)
  {
    this.variables = new LinkedHashSet<>(Arrays.asList(variables));
  }
  
  public LocAtCollisionsProcessor(Collection<RealVariable> variables)
  {
    this.variables = new LinkedHashSet<RealVariable>(variables);
  }

  @Override
  public void init(LocalRFSampler sampler)
  {
    this.sampler = sampler;
    update(0.0);
  }

  @Override
  public void processRay(RealVariable var, TrajectoryRay ray, double time,
      LocalRFSampler sampler)
  {
    if (time == prevTime)
      return;
    
    if (!variables.contains(var))
      return;
    
    update(time);
    
    prevTime = time;
  }
  
  private void update(double time) 
  {
    DoubleMatrix current = new DoubleMatrix(variables.size());
    int i = 0;
    for (RealVariable aVar : variables)
    {
      sampler.updateVariable(aVar, time);
      current.put(i, aVar.getValue());
      i++;
    }
    locationsAtCollisions.add(current);
  }
  
}