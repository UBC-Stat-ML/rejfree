package rejfree;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;

import com.google.common.collect.Lists;

import bayonet.distributions.Exponential;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.LBFGSMinimizer;
import briefj.opt.Option;

import static rejfree.StaticUtils.*;

public class SimpleRFSampler
{
  /**
   * *Negative* log density of the target distribution.
   */
  private final DifferentiableFunction energy;
  private final PegasusConvexCollisionSolver solver = new PegasusConvexCollisionSolver();
  private final SimpleRFSamplerOptions options;
  
  private List<DoubleMatrix> trajectory = Lists.newArrayList();
  private List<DoubleMatrix> samples = Lists.newArrayList();
  
  private DoubleMatrix currentPosition, currentVelocity;

  public static class SimpleRFSamplerOptions
  {
    @Option
    public double refreshRate = 0.0001;
    
    @Option
    public double collectRate = 1.0;

    @Option
    public boolean useInformedVelocityUpdate = false;
  }
  
  private SummaryStatistics collisionToRefreshmentRatio = new SummaryStatistics();
  private SummaryStatistics collectedPerEvent = new SummaryStatistics();
  
  /**
   * 
   * @param energy The negative log density of the target, assumed to be convex
   */
  public SimpleRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition, SimpleRFSamplerOptions options)
  {
    this.energy = energy;
    this.options = options;
    this.currentPosition = initialPosition;
    this.currentVelocity = uniformOnUnitBall(energy.dimension(), new Random(1));
  }
  
  public SimpleRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition)
  {
    this(energy, initialPosition, new SimpleRFSamplerOptions());
  }
  
  public static SimpleRFSampler initializeRFWithLBFGS(DifferentiableFunction energy, SimpleRFSamplerOptions options)
  {
    return new SimpleRFSampler(energy, optimizePosition(energy), options);
  }
  
  public static SimpleRFSampler initializeRFWithLBFGS(DifferentiableFunction energy)
  {
    return initializeRFWithLBFGS(energy, new SimpleRFSamplerOptions());
  }
  
  /**
   * 
   * @param dimension
   * @param rand
   * @return A random vector of unit norm.
   */
  public static DoubleMatrix uniformOnUnitBall(int dimension, Random rand)
  {
    DoubleMatrix random = new DoubleMatrix(dimension);
    for (int i = 0; i < dimension; i++)
      random.data[i] = rand.nextGaussian();
    double norm = random.norm2();
    return random.muli(1.0/norm);
  }

  private static DoubleMatrix optimizePosition(DifferentiableFunction energy)
  {
    double [] min = new LBFGSMinimizer().minimize(energy, new DoubleMatrix(energy.dimension()).data, 1e-2);
    return new DoubleMatrix(min);
  }
  private DoubleMatrix _cachedOptimizePosition = null;
  private DoubleMatrix cachedOptimizePosition()
  {
    if (_cachedOptimizePosition == null)
      _cachedOptimizePosition = optimizePosition(this.energy);
    return _cachedOptimizePosition;
  }

  public void iterate(Random rand, int numberOfIterations)
  {
    trajectory.add(currentPosition);
    
    for (int iter = 0; iter < numberOfIterations; iter++)
    {
      // simulate event
      final double exponential = - Math.log(rand.nextDouble());
      double collisionTime = solver.collisionTime(currentPosition, currentVelocity, energy, exponential);
      double refreshTime = Exponential.generate(rand, options.refreshRate);
      double eventTime = Math.min(collisionTime, refreshTime);
      collisionToRefreshmentRatio.addValue(collisionTime/refreshTime);
      
      // collect state
      collectSamples(currentPosition, currentVelocity, eventTime, rand);
      
      // update state
      boolean collisionOccurs = collisionTime < refreshTime;
      currentPosition = position(currentPosition, currentVelocity, eventTime);
      trajectory.add(currentPosition);
      if (collisionOccurs)
        currentVelocity = Bouncer.bounce(currentVelocity, gradient(currentPosition));
      else
        currentVelocity = refreshVelocity(currentPosition, currentVelocity, rand); 
    }
  }
  

  private DoubleMatrix refreshVelocity(DoubleMatrix currentPosition,
      DoubleMatrix currentVelocity, Random rand)
  {
    if (rand.nextBoolean() || !options.useInformedVelocityUpdate)
      return uniformOnUnitBall(currentVelocity.length, rand);
    else
    {
      DoubleMatrix difference = currentPosition.sub(cachedOptimizePosition());
      difference.muli(1.0/difference.norm2());
      if (rand.nextBoolean())
        difference.muli((-1));
      return difference;
    }
  }

  private void collectSamples(DoubleMatrix initialPosition,
      DoubleMatrix velocity, double eventTime, Random rand)
  {
    double timeConsumed = Exponential.generate(rand, options.collectRate);
    int nCollected = 0;
    while (timeConsumed < eventTime)
    {
      nCollected++;
      samples.add(position(initialPosition, velocity, timeConsumed));
      timeConsumed += Exponential.generate(rand, options.collectRate);
    }
    collectedPerEvent.addValue(nCollected);
  }

  private DoubleMatrix gradient(DoubleMatrix position)
  {
    return new DoubleMatrix(energy.derivativeAt(position.data));
  }
  
  public List<DoubleMatrix> getTrajectory()
  {
    return trajectory;
  }

  public SummaryStatistics getCollisionToRefreshmentRatio()
  {
    return collisionToRefreshmentRatio;
  }

  public List<DoubleMatrix> getSamples()
  {
    return samples;
  }

  public SummaryStatistics getCollectedPerEvent()
  {
    return collectedPerEvent;
  }

  public int dimensionality()
  {
    return energy.dimension();
  }

  public DoubleMatrix getCurrentPosition()
  {
    return currentPosition;
  }

  public void setCurrentPosition(DoubleMatrix currentPosition)
  {
    this.currentPosition = currentPosition;
  }
}
