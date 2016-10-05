package rejfree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import bayonet.distributions.Poisson;

import static java.lang.Math.*;

public class TestAdaptiveThinning
{
  static interface PP
  {
    double intensity(double t);
  }
  
  static interface ProposalPP
  {
    Optional<Double> sample(Random rand, double s);
    double delta(double s);
    double intensity(double s, double t);
  }
  
  static double adaptiveThinningSampler(Random rand, ProposalPP proposal, PP target)
  {
    double s = 0.0;
    do
    {
      double delta = proposal.delta(s);
      double tau = proposal.sample(rand, s).orElse(s + delta);
      if (s + delta <= tau)
        s = s + delta;
      else
      {
        if (rand.nextDouble() > target.intensity(tau) / proposal.intensity(s, tau)) 
          s = tau; 
        else 
          return tau;
      }
    }
    while (true);
  }
  
  // Example
  
  static class SineTarget implements PP
  {
    @Override
    public double intensity(double s)
    {
      return sin(s) + 1.0;
    }
  }
  
  /********
   
   +-+
   | +-+
   +---+
   s   s+2
   
   *********/
  static class TetrisProposal implements ProposalPP
  {

    @Override
    public double intensity(double s, double t)
    {
           if (t < s)       return 0.0;
      else if (t > s + 2.0) return 0.0;
      else if (t < s + 1.0) return 2.0;
      else                  return 1.0;
    }

    @Override
    public Optional<Double> sample(Random rand, double s)
    {
      int nPoints = Poisson.generate(rand, 3);
      List<Double> points = new ArrayList<>();
      for (int i = 0; i < nPoints; i++)
      {
        int section = rand.nextInt(3);
        double point = rand.nextDouble();
        points.add(s + point + (section <= 1 ? 0.0 : 1.0));
      }
      Collections.sort(points);
      if (points.isEmpty())
        return Optional.empty();
      else
        return Optional.of(points.get(0));
    }

    @Override
    public double delta(double s)
    {
      return s + 1.0;
    }
  }
  
  @Test
  public void test()
  {
    Random rand = new Random(1);
    int nSamples = 1_000_000;
    
    double precision = 1/sqrt(nSamples) * 10;
    
    List<Double> samples = new ArrayList<>();
    ProposalPP proposal = new TetrisProposal();
    PP target = new SineTarget();
    for (int i = 0; i < nSamples; i++)
      samples.add(adaptiveThinningSampler(rand, proposal, target));
    
    // check the CDF against analytical for a few grid points
    for (double grid = 0.0; grid < 10; grid += 0.2)
    {
      final double currentGrid = grid;
      double empirical = ((double) samples.stream().filter(sample -> sample < currentGrid).count()) / (double) nSamples;
      double analytic = 1.0 - exp(cos(currentGrid) - currentGrid - 1.0);
      Assert.assertEquals(analytic, empirical, precision);
    }
  }
}
