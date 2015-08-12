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

public class GlobalRFSampler
{
  /**
   * *Negative* log density of the target distribution.
   */
  private final DifferentiableFunction energy;
  private final PegasusConvexCollisionSolver solver = new PegasusConvexCollisionSolver();
  private final RFSamplerOptions options;
  
  private List<DoubleMatrix> trajectory = Lists.newArrayList();
  private List<DoubleMatrix> samples = Lists.newArrayList();
  
  private DoubleMatrix currentPosition, currentVelocity;

  public static class RFSamplerOptions
  {
    @Option
    public double refreshRate = 1.0;
    
    @Option
    public double collectRate = 1.0;

    @Option
    public boolean useInformedVelocityUpdate = false;

    @Option
    public boolean useLocalRefreshment = false;
  }
  
  private SummaryStatistics collisionToRefreshmentRatio = new SummaryStatistics();
  private SummaryStatistics collectedPerEvent = new SummaryStatistics();
  
  /**
   * 
   * @param energy The negative log density of the target, assumed to be convex
   */
  public GlobalRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition, RFSamplerOptions options)
  {
    this.energy = energy;
    this.options = options;
    this.currentPosition = initialPosition;
    this.currentVelocity = null; 
  }
  
  public GlobalRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition)
  {
    this(energy, initialPosition, new RFSamplerOptions());
  }
  
  public static GlobalRFSampler initializeRFWithLBFGS(DifferentiableFunction energy, RFSamplerOptions options)
  {
    return new GlobalRFSampler(energy, optimizePosition(energy), options);
  }
  
  public static GlobalRFSampler initializeRFWithLBFGS(DifferentiableFunction energy)
  {
    return initializeRFWithLBFGS(energy, new RFSamplerOptions());
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
  
  private DoubleMatrix 
    mean, 
    variance;
  
  public DoubleMatrix getMean()
  {
    return mean;
  }
  
  public DoubleMatrix getVariance()
  {
    return variance;
  }

  public void iterate(Random rand, int numberOfIterations)
  {
    if (currentVelocity == null)
      currentVelocity = uniformOnUnitBall(energy.dimension(), rand);
    
    trajectory.add(currentPosition);
    double totalTime = 0.0;
    mean = new DoubleMatrix(dimensionality());
    variance = new DoubleMatrix(dimensionality(),dimensionality());
    for (int iter = 0; iter < numberOfIterations; iter++)
    {
      // simulate event
      final double exponential = StaticUtils.generateUnitRateExponential(rand);
      double collisionTime = solver.collisionTime(currentPosition, currentVelocity, energy, exponential);
      double refreshTime = options.refreshRate == 0 ? Double.POSITIVE_INFINITY : Exponential.generate(rand, options.refreshRate);
      double eventTime = Math.min(collisionTime, refreshTime);
      totalTime += eventTime;
      collisionToRefreshmentRatio.addValue(collisionTime/refreshTime);
      
      // collect state
      collectSamples(currentPosition, currentVelocity, eventTime, rand);
      
      // update state
      boolean collisionOccurs = collisionTime < refreshTime;
      currentPosition = position(currentPosition, currentVelocity, eventTime);
      trajectory.add(currentPosition);
      if (collisionOccurs)
        currentVelocity = StaticUtils.bounce(currentVelocity, gradient(currentPosition));
      else
        currentVelocity = refreshVelocity(currentPosition, currentVelocity, rand); 
    }
    mean.divi(totalTime);
    variance.divi(totalTime);
  }
  
  public void setVelocity(DoubleMatrix velocity)
  {
    this.currentVelocity = velocity.dup();
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
    mean.addi(initialPosition.mul(eventTime)
        .addi(velocity .mul(eventTime*eventTime/2.0)));
    variance
        .addi(initialPosition.mmul(initialPosition.transpose()).mul(eventTime))
        .addi(initialPosition.mmul(velocity.transpose()).mul(eventTime*eventTime))
        .addi(velocity.mmul(velocity.transpose()).mul(eventTime*eventTime*eventTime/3.0));
    
    if (options.collectRate == 0.0)
      return;
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
