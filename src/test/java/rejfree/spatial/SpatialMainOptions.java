package rejfree.spatial;

import java.io.File;
import java.util.Random;

import briefj.opt.Option;



public class SpatialMainOptions
{
  @Option
  public int nSamples = 1000;
  
  @Option
  public double drift = 0.2;
  
  @Option
  public double init = 10;
  
  @Option
  public Random random = new Random(1);
  
  @Option(required = true)
  public File dataFolder;
  
  public File getRDataFile()
  {
    if (init != 10 || drift != 0.2) // TODO: right now, hard-coded in stan code
      throw new RuntimeException();
    
    return new File(dataFolder, PreprocessData.R_OUTPUT_NAME);
  }
  
  public File getGeoDataCSVFile()
  {
    return new File(dataFolder, PreprocessData.GEO_OUTPUT_NAME);
  }
  
  public File getAccidentsDataCSVFile()
  {
    return new File(dataFolder, PreprocessData.ACCIDENTS_OUTPUT_NAME);
  }
}
