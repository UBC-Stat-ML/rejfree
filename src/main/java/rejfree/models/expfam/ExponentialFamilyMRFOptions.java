package rejfree.models.expfam;


import briefj.opt.Option;



public class ExponentialFamilyMRFOptions
{
  @Option
  public double diag = 1;
  
  @Option
  public double offDiag = 0.5;
  
  @Option
  public int nRows = 10;
  
  @Option
  public int nCols = 10;
  
  public boolean hasLikelihood()
  {
    return true;
  }
}