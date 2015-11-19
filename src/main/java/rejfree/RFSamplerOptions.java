package rejfree;

import briefj.opt.Option;



public class RFSamplerOptions
{
  @Option
  public double refreshRate = 1.0;
  
  @Option
  public double collectRate = 1.0;

  @Option
  public boolean useInformedVelocityUpdate = false;

  @Option
  public boolean useLocalRefreshment = true;
  
  @Option
  public boolean restrictVelocityNorm = true;

  @Option
  public boolean usePartialRefreshment = false;

  @Option
  public double alpha = 1.0;
  
  @Option
  public double beta = 4.0;
}