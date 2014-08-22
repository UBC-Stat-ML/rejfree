package rejfree;

import java.util.Random;





public class NormalExample
{
  
  public static void main(String [] args)
  {
    NormalEnergy energy = NormalEnergy.isotropic(2);
    SimpleRFSampler sampler = new SimpleRFSampler(energy);
    Random rand = new Random(1);
    sampler.iterate(rand, 100);
  }
}
