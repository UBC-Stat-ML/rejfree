package rejfree;

import briefj.opt.Option;



public class RFSamplerOptions
{
  @Option
  public double refreshRate = 1.0;
  
  @Option
  public double collectRate = 1.0;
  
  @Option
  public RefreshmentMethod refreshmentMethod = RefreshmentMethod.LOCAL;
  
  public static enum RefreshmentMethod
  {
    GLOBAL, LOCAL, RESTRICTED, PARTIAL;
  }

  @Option
  public double alpha = 1.0;
  
  @Option
  public double beta = 4.0;
}