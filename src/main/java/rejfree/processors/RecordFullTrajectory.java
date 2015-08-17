package rejfree.processors;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;

import rejfree.local.LocalRFSampler;
import rejfree.local.TrajectoryRay;
import blang.variables.RealVariable;
import briefj.Indexer;



public class RecordFullTrajectory implements RayProcessor
{
  public final List<DoubleMatrix> samples = new ArrayList<>();
  public final Indexer<RealVariable> variablesIndexer;
  private double lastT = 0.0;
  
  public RecordFullTrajectory(Indexer<RealVariable> variablesIndexer)
  {
    super();
    this.variablesIndexer = variablesIndexer;
  }

  @Override
  public void processRay(RealVariable var, TrajectoryRay ray,
      double time, LocalRFSampler sampler)
  {
    if (time == lastT)
      return;
    sampler.updateAllVariables(time);
    process(sampler);
    lastT = time;
  }

  @Override
  public void init(LocalRFSampler sampler)
  {
    process(sampler);
  }
  
  private void process(LocalRFSampler sampler)
  {
    DoubleMatrix current = new DoubleMatrix(variablesIndexer.size());
    int i = 0;
    for (RealVariable var : variablesIndexer.objectsList())
      current.put(i++, var.getValue());
    samples.add(current);
  }
}