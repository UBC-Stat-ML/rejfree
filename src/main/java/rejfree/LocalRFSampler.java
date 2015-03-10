package rejfree;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;

import org.jblas.DoubleMatrix;

import rejfree.SimpleRFSampler.SimpleRFSamplerOptions;
import bayonet.distributions.Exponential;
import blang.ProbabilityModel;
import blang.factors.Factor;
import blang.variables.RealVariable;



public class LocalRFSampler
{
  private final EventQueue<CollisionFactor> collisionQueue = new EventQueue<>();
  private final Map<RealVariable, TrajectoryRay> trajectories = new HashMap<>();
  
  private final ProbabilityModel model;
  private final SimpleRFSamplerOptions simpleRFOptions;
  
  public LocalRFSampler(ProbabilityModel model, SimpleRFSamplerOptions options)
  {
    this.model = model;
    this.simpleRFOptions = options;
    globalVelocityRefreshment(new Random(1), 0.0);
  }
  
  private void globalVelocityRefreshment(Random rand, double refreshmentTime)
  {
    final List<Object> variables = model.getLatentVariables();
    final int dimensionality = variables.size();
    final DoubleMatrix uniformOnUnitBall = SimpleRFSampler.uniformOnUnitBall(variables.size(), rand);
    for (int i = 0; i < dimensionality; i++)
    {
      RealVariable variable = (RealVariable) variables.get(i);
      double currentVelocity = uniformOnUnitBall.get(i);
      if (trajectories.containsKey(variable))
        updateTrajectory(refreshmentTime, variable, currentVelocity);
      else
        trajectories.put(variable, new TrajectoryRay(refreshmentTime, 0.0, uniformOnUnitBall.get(i)));
    }
    for (Factor factor : model.linearizedFactors())
      updateCandidateCollision(rand, (CollisionFactor) factor);
  }

  public void iterate(Random rand, int numberOfIterations)
  {
    double nextGlobalRefreshmentTime = Exponential.generate(rand, simpleRFOptions.refreshRate);
    for (int iter = 0; iter < numberOfIterations; iter++)
    {
      double nextCollisionTime = collisionQueue.peekTime();
      if (nextCollisionTime < nextGlobalRefreshmentTime)
        doCollision(rand);
      else
      {
        globalVelocityRefreshment(rand, nextGlobalRefreshmentTime);
        nextGlobalRefreshmentTime += Exponential.generate(rand, simpleRFOptions.refreshRate);
      }
    }
  }
  
  /**
   * Perform one collision and the associated updates to the queue
   * @param rand
   */
  private <CollisionType> void doCollision(Random rand)
  {
    // 0 - pop a collision factor
    final Entry<Double,CollisionFactor> collision = collisionQueue.pollEvent();
    final double collisionTime = collision.getKey();
    final CollisionFactor collisionFactor = collision.getValue();
    
    Collection<?> immediateNeighborVariables = model.neighborVariables(collisionFactor);
    Collection<CollisionFactor> neighborFactors = neighborFactors(immediateNeighborVariables);
    Collection<?> extendedNeighborVariables = neighborVariables(neighborFactors);
    
    // 1- update RealVariables in extended neighborhood
    for (Object variable : extendedNeighborVariables)
      updateVariable(variable, collisionTime);
    
    // 2- update rays for variables in immediate neighborhood (and process)
    updateTrajectories(collisionFactor, collisionTime);
    
    // 3- recompute the collisions for the other factors touching the variables (including the one we just popped)
    for (CollisionFactor factor : neighborFactors)
      updateCandidateCollision(rand, factor);
  }
  
  private void updateCandidateCollision(Random rand, CollisionFactor factor)
  {
    collisionQueue.remove(factor);
    double exponentialRealization = Exponential.generate(rand, 1.0);
    double candidateCollisionTime = factor.getCollisionTime(exponentialRealization);
    collisionQueue.add(factor, candidateCollisionTime);
  }
  
  /**
   * @param neighborFactors
   * @return Distinct neighbors variables of the provided factors
   */
  private Collection<?> neighborVariables(
      Collection<CollisionFactor> neighborFactors)
  {
    HashSet<Object> result = new LinkedHashSet<>();
    for (CollisionFactor factor : neighborFactors)
      for (Object variable : model.neighborVariables(factor))
        result.add(variable);
    return result;
  }

  /**
   * @param immediateNeighborVariables
   * @return Distinct neighbor factors of the provided variables
   */
  private Collection<CollisionFactor> neighborFactors(
      Collection<?> immediateNeighborVariables)
  {
    HashSet<CollisionFactor> result = new LinkedHashSet<>();
    for (Object variable : immediateNeighborVariables)
      for (Factor factor : model.neighborFactors(variable))
        result.add((CollisionFactor) factor);
    return result;
  }

  /**
   * Update all trajectories affected by one collision (i.e. those connected to the
   * colliding factor)
   * @param collisionFactor
   * @param collisionTime
   */
  private void updateTrajectories(CollisionFactor collisionFactor, double collisionTime)
  {
    DoubleMatrix gradient = collisionFactor.gradient();
    DoubleMatrix oldVelocity = getVelocityMatrix(collisionFactor, gradient.length);
    DoubleMatrix newVelocity = Bouncer.bounce(oldVelocity, gradient);
    
    final int length = newVelocity.length;
    for (int i = 0; i < length; i++)
    {
      RealVariable variable = collisionFactor.getVariable(i);
      double newVelocityCoordinate = newVelocity.get(i);
      updateTrajectory(collisionTime, variable, newVelocityCoordinate);
    }
  }

  /**
   * Update one trajectory coordinate (a single real variable)
   * @param time
   * @param variable
   * @param newVelocity
   */
  private void updateTrajectory(double time, RealVariable variable,
      double newVelocity)
  {
    TrajectoryRay oldRay = trajectories.get(variable);
    double newPosition = oldRay.position(time);
    TrajectoryRay newRay = new TrajectoryRay(time, newPosition, newVelocity);
    trajectories.put(variable, newRay);
    processRay(oldRay, time);
  }

  private void processRay(TrajectoryRay ray, double timeTheRayEnds)
  {
    // TODO Auto-generated method stub
  }

  private DoubleMatrix getVelocityMatrix(CollisionFactor factor, int length)
  {
    DoubleMatrix result = new DoubleMatrix(length);
    for (int i = 0; i < length; i++)
      result.put(i, trajectories.get(factor.getVariable(i)).velocity_t);
    return result;
  }

  private void updateVariable(Object _variable, double currentTime)
  {
    RealVariable variable = (RealVariable) _variable;
    TrajectoryRay ray = trajectories.get(variable);
    variable.setValue(ray.position(currentTime));
  }
}
