package rejfree.local;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import com.google.common.base.Stopwatch;

import rejfree.StaticUtils;
import rejfree.GlobalRFSampler.RFSamplerOptions;
import bayonet.distributions.Exponential;
import blang.ProbabilityModel;
import blang.MCMCFactory.Factories;
import blang.MCMCFactory.MCMCOptions;
import blang.factors.Factor;
import blang.processing.NodeProcessorFactory;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.processing.ProcessorFactory;
import blang.variables.RealVariable;
import briefj.Indexer;
import briefj.collections.Counter;



public class LocalRFSampler
{
  private final EventQueue<CollisionFactor> _collisionQueue = new EventQueue<>();
  private final Map<CollisionFactor,Boolean> isCollisionMap = new LinkedHashMap<>();
  private final Map<RealVariable, TrajectoryRay> trajectories = new HashMap<>();
  private final List<CollisionFactor> allFactors;
  
  public final ProbabilityModel model;
  private final RFSamplerOptions rfOptions;
  public final MCMCOptions mcmcOptions;
  
  public final List<Processor> processors = new ArrayList<Processor>();
  public final List<RayProcessor> rayProcessors = new ArrayList<>();
  
  private int nCollisions = 0;
  private int nCollidedVariables = 0;
  
  private int nRefreshments = 0;
  private int nRefreshedVariables = 0;
  
  public LocalRFSampler(ProbabilityModel model, RFSamplerOptions options)
  {
    this.model = model;
    this.rfOptions = options;
    this.mcmcOptions = new MCMCOptions();
    mcmcOptions.burnIn = 0;
    mcmcOptions.thinningPeriod = 1;
    mcmcOptions.nMCMCSweeps = Integer.MAX_VALUE;
    // mcmcOptions.progressCODA = true;  <-- avoid this, it makes things slow
    
    allFactors = new ArrayList<>();
    for (Factor f : model.linearizedFactors())
      allFactors.add((CollisionFactor) f);
  }
  
  public void addPointProcessor(Processor processor)
  {
    this.processors.add(processor);
  }
  
  public void addDefaultPointProcessors()
  {
    Factories<ProcessorFactory,NodeProcessorFactory> processorFactories = new Factories<ProcessorFactory,NodeProcessorFactory>(new NodeProcessorFactory());
    for (ProcessorFactory f : processorFactories.factories)
      processors.addAll(f.build(model));
  }
  
  public MomentRayProcessor addDefaultMomentRayProcessor()
  {
    MomentRayProcessor momentProcessor = new MomentRayProcessor();
    rayProcessors.add(momentProcessor);
    return momentProcessor;
  }
  
  private void processRay(RealVariable var, TrajectoryRay ray, double timeTheRayEnds)
  {
    for (RayProcessor processor : rayProcessors)
      processor.processRay(var, ray, timeTheRayEnds, this);
  }
  
  @SuppressWarnings({"rawtypes", "unchecked"})
  public RecordFullTrajectory addRecordFullTrajectoryProcessor()
  {
    Indexer variablesIndexer = new Indexer<>(model.getLatentVariables());
    RecordFullTrajectory processor = new RecordFullTrajectory(variablesIndexer);
    rayProcessors.add(processor);
    return processor;
  }
  
  private int pointCollectIter = 0;
  private void processPoint()
  {
    for (Processor p : processors)
      p.process(new ProcessorContext(pointCollectIter++, model, mcmcOptions));
  }
  
  private void globalVelocityRefreshment(Random rand, double refreshmentTime, boolean initializing)
  {
    final List<Object> variables = model.getLatentVariables();
    final int dimensionality = variables.size();
    final DoubleMatrix uniformOnUnitBall = StaticUtils.uniformOnUnitBall(variables.size(), rand);
    
    nRefreshments++;
    nRefreshedVariables += variables.size();
    
    for (int i = 0; i < dimensionality; i++)
    {
      RealVariable variable = (RealVariable) variables.get(i);
      double currentVelocity = uniformOnUnitBall.get(i);
      if (initializing)
      {
        initTrajectory(refreshmentTime, variable, currentVelocity);
        updateVariable(variable, refreshmentTime);
      }
      else
      {
        updateVariable(variable, refreshmentTime);
        updateTrajectory(refreshmentTime, variable, currentVelocity);
      }
    }
    
    for (Factor factor : model.linearizedFactors())
      updateCandidateCollision(rand, (CollisionFactor) factor, refreshmentTime);
  }
  
  private void localVelocityRefreshment(Random rand, double refreshmentTime)
  {
    // sample a factor
    CollisionFactor f = allFactors.get(rand.nextInt(allFactors.size()));
    Collection<?> immediateNeighborVariables = model.neighborVariables(f);
    Collection<CollisionFactor> neighborFactors = neighborFactors(immediateNeighborVariables);
    Collection<?> extendedNeighborVariables = neighborVariables(neighborFactors);
    
    nRefreshments++;
    nRefreshedVariables += immediateNeighborVariables.size();
    
    // compute velocity squared norm
    double sum = 0.0;
    for (Object variable : immediateNeighborVariables)
    {
      double currentVComponent = trajectories.get(variable).velocity_t;
      sum += currentVComponent * currentVComponent;
    }
    
    // sample new velocity vector, rescale it
    final DoubleMatrix newVelocity = StaticUtils.uniformOnUnitBall(immediateNeighborVariables.size(), rand);
    newVelocity.muli(Math.sqrt(sum));
    
    for (Object variable : extendedNeighborVariables)
      updateVariable(variable, refreshmentTime);
    
    // 2- update rays for variables in immediate neighborhood (and process)
    int d = 0;
    for (Object variable : immediateNeighborVariables)
      updateTrajectory(refreshmentTime, (RealVariable) variable, newVelocity.get(d++));
    
    // 3- recompute the collisions for the other factors touching the variables (including the one we just popped)
    for (CollisionFactor factor : neighborFactors)
      updateCandidateCollision(rand, factor, refreshmentTime);
  }

  private void initTrajectory(double refreshmentTime, RealVariable variable, double currentVelocity)
  {
    if (trajectories.containsKey(variable))
      throw new RuntimeException();
    trajectories.put(variable, new TrajectoryRay(refreshmentTime, variable.getValue(), currentVelocity));
  }
  
  public void iterate(Random rand, int maxNumberOfIterations)
  {
    iterate(rand, maxNumberOfIterations, Double.POSITIVE_INFINITY, Long.MAX_VALUE);
  }
  
  public void iterate(Random rand, int maxNumberOfIterations, double maxTrajectoryLen)
  {
    iterate(rand, maxNumberOfIterations, maxTrajectoryLen, Long.MAX_VALUE);
  }

  private double currentTime = 0.0;
  public void iterate(Random rand, int maxNumberOfIterations, double maxTrajectoryLen, long maxTimeMilli)
  {
    Stopwatch watch = maxTimeMilli == Long.MAX_VALUE ? null : Stopwatch.createStarted();
    if (currentTime == 0.0)
    {
      globalVelocityRefreshment(rand, 0.0, true);
      for (RayProcessor rayProc : rayProcessors)
        rayProc.init(this);
    }
    double nextRefreshmentTime = rfOptions.refreshRate == 0 ? Double.POSITIVE_INFINITY : Exponential.generate(rand, rfOptions.refreshRate);
    mainLoop : for (int iter = 0; iter < maxNumberOfIterations; iter++)
    {
      if (watch != null && watch.elapsed(TimeUnit.MILLISECONDS) > maxTimeMilli)
        break mainLoop;
      
      double nextCollisionTime = _collisionQueue.peekTime(); 
      double nextEventTime = Math.min(nextCollisionTime, nextRefreshmentTime);
      if (nextEventTime > maxTrajectoryLen)
      {
        currentTime = maxTrajectoryLen;
        break mainLoop;
      }
      collectSamples(currentTime, nextEventTime, rand);
      if (nextCollisionTime < nextRefreshmentTime)
      {
        doCollision(rand);
        currentTime = nextCollisionTime;
      }
      else
      {
        currentTime = nextRefreshmentTime;
        if (rfOptions.useLocalRefreshment)
          localVelocityRefreshment(rand, nextRefreshmentTime);
        else
          globalVelocityRefreshment(rand, nextRefreshmentTime, false);
        nextRefreshmentTime += Exponential.generate(rand, rfOptions.refreshRate);
      }
    }
    // Ensure that the state is globally consistent at the end
    updateAllVariables(currentTime);
    
    // Ensure that the remaining rays are processed 
    globalVelocityRefreshment(rand, currentTime, false); 
  }
  
  private void collectSamples(double currentTime, double nextEventTime, Random rand)
  {
    if (rfOptions.collectRate == 0.0)
      return;
    double timeConsumed = currentTime + Exponential.generate(rand, rfOptions.collectRate);
    while (timeConsumed < nextEventTime)
    {
      updateAllVariables(timeConsumed);
      processPoint();
      timeConsumed += Exponential.generate(rand, rfOptions.collectRate);
    }
  }
  
  /**
   * Perform one collision and the associated updates to the queue
   * @param rand
   */
  private <CollisionType> void doCollision(Random rand)
  {
    // 0 - pop a collision factor
    final Entry<Double,CollisionFactor> collision = _collisionQueue.pollEvent();
    
    final double collisionTime = collision.getKey();
    final CollisionFactor collisionFactor = collision.getValue();
    final boolean isActualCollision = isCollisionMap.get(collisionFactor);
    
    // TODO: if isActualCollision is false, it might be possible to update only a subset
    // of the nodes below
    Collection<?> immediateNeighborVariables = model.neighborVariables(collisionFactor);
    Collection<CollisionFactor> neighborFactors = neighborFactors(immediateNeighborVariables);
    Collection<?> extendedNeighborVariables = neighborVariables(neighborFactors);
    
    nCollisions++;
    nCollidedVariables += immediateNeighborVariables.size();
    
    // 1- update RealVariables in extended neighborhood
    for (Object variable : extendedNeighborVariables)
      updateVariable(variable, collisionTime);
    
    // 2- update rays for variables in immediate neighborhood (and process)
    if (isActualCollision)
      collideTrajectories(collisionFactor, collisionTime);
    
    if (isActualCollision)
    {
      // 3- recompute the collisions for the other factors touching the variables (including the one we just popped)
      for (CollisionFactor factor : neighborFactors)
        updateCandidateCollision(rand, factor, collisionTime);
    }
    else
    {
      // 3b- the collision is actually just a trigger to recompute the next collision time
      updateCandidateCollision(rand, collisionFactor, collisionTime);
    }
  }
  
  private void updateCandidateCollision(Random rand, CollisionFactor factor, double currentTime)
  {
    _collisionQueue.remove(factor);
    
    CollisionContext context = new CollisionContext(rand, getVelocityMatrix(factor));
    Pair<Double, Boolean> collisionInfo = factor.getLowerBoundForCollisionDeltaTime(context);
    
    double candidateCollisionTime = currentTime + collisionInfo.getLeft();
    isCollisionMap.put(factor, collisionInfo.getRight());
    
    if (_collisionQueue.containsTime(candidateCollisionTime))
    {
      System.err.println("The sampler has hit an event of probability zero: two collisions scheduled exactly at the same time.");
      System.err.println("Because of numerical precision, this could possibly happen, but very rarely.");
      
      System.err.println("For internal implementation reasons, one of the collisions at time " + candidateCollisionTime + " was moved to " + (candidateCollisionTime + epsilon));
      candidateCollisionTime += epsilon;
    }
    
    _collisionQueue.add(factor, candidateCollisionTime);
  }
  
  private static final double epsilon = 1e-10;
  
  /**
   * @param neighborFactors
   * @return Distinct neighbors variables of the provided factors
   */
  private Collection<?> neighborVariables(
      Collection<CollisionFactor> neighborFactors)
  {
    HashSet<Object> result = new LinkedHashSet<>();
    for (CollisionFactor factor : neighborFactors)
      for (Object variable : model.neighborLatentVariables(factor))
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
  private void collideTrajectories(CollisionFactor collisionFactor, double collisionTime)
  {
    DoubleMatrix gradient = collisionFactor.gradient();
    DoubleMatrix oldVelocity = getVelocityMatrix(collisionFactor);
    DoubleMatrix newVelocity = StaticUtils.bounce(oldVelocity, gradient);
    
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
    processRay(variable, oldRay, time);
  }

  private DoubleMatrix getVelocityMatrix(CollisionFactor factor)
  {
    final int length = factor.nVariables();
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
  
  public void updateAllVariables(double currentTime)
  {
    for (Object var : model.getLatentVariables())
      updateVariable(var, currentTime);
  }
  
  public int getNCollisions()
  {
    return nCollisions;
  }
  
  public int getNCollidedVariables()
  {
    return nCollidedVariables;
  }
  
  public int getNRefreshments()
  {
    return nRefreshments;
  }
  
  public int getNRefreshedVariables()
  {
    return nRefreshedVariables;
  }

  public static class MomentRayProcessor implements RayProcessor
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
    
    private double indefIntegralForMean(double x0, double v, double t)
    {
      return x0 * t + v * t*t / 2.0;
    }
    
    private double indefIntegralForVar(double x0, double v, double t)
    {
      return x0*x0 * t + x0 * v * t*t + v*v * t*t*t / 3.0;
    }
  };
}
