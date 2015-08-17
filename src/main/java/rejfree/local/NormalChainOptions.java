package rejfree.local;

import java.util.Random;

import briefj.opt.Option;



public class NormalChainOptions
{
  @Option
  public double diag = 1;
  
  @Option
  public double offDiag = 0.5;
  
  @Option
  public int nPairs = 3; 
  
  @Option
  public Random random = new Random(1);

  @Option
  public boolean useAnalytic = true;
  
  @Option(gloss = "Use the local structure to create several sparse normal precision factors")
  public boolean useLocal = true;
}