package rejfree.spatial;

import java.io.File;

import bayonet.rplot.RJavaBridge;
import bayonet.rplot.RUtils;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class PostprocessIntersectionSamples extends RJavaBridge implements Runnable
{
  @Option(required = true)
  public File samples;
  
  public File getIntensitiesOutput()
  {
    return Results.getFileInResultFolder("intensities.pdf");
  }
  
  private PostprocessIntersectionSamples() {}
  
  public PostprocessIntersectionSamples(File samples)
  {
    this.samples = samples;
  }
  
//  @Option(required = true)
//  public File streetAtCornerCsvFile;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new PostprocessIntersectionSamples());
  }
  
  @Override
  public void run()
  {
    RUtils.callRBridge(this);
  }
//
//  @Override
//  public void run()
//  {
//    Map<Integer,String>  index = readIndex(streetAtCornerCsvFile);
//    Map<Integer,SummaryStatistics> statistics = new LinkedHashMap<>();
//    for (Map<String,String> line : BriefIO.readLines(samples).indexCSV('#'))
//    {
//      for (Integer streetAtCornerId : index.keySet())
//      {
//        String stanName = "logIntensity." + (streetAtCornerId+1);
//        double sample = Math.exp(Double.parseDouble(line.get(stanName)));
//        BriefMaps.getOrPut(statistics, streetAtCornerId, new SummaryStatistics()).addValue(sample);
//      }
//    }
//    for (Integer streetAtCornerId : index.keySet())
//    {
//      String readableName = index.get(streetAtCornerId);
//      SummaryStatistics stats = statistics.get(streetAtCornerId);
//      System.out.println(readableName + "\t" + stats.getMean() + "\t" + stats.getStandardDeviation());
//    }
//  }
//
//
//  private static Map<Integer, String> readIndex(File streetAtCornerCsvFile2)
//  {
//    Map<Integer,String> result = new LinkedHashMap<Integer, String>();
//    
//    for (Map<String,String> line : BriefIO.readLines(streetAtCornerCsvFile2).indexCSV())
//    {
//      addToIndex(result, line.get("streetAtCorner1"), line.get("streetAtCorner1Index"));
//      addToIndex(result, line.get("streetAtCorner2"), line.get("streetAtCorner2Index"));
//    }
//    return result;
//  }
//
//  private static void addToIndex(
//      Map<Integer, String> result, 
//      String name,
//      String indexStr)
//  {
//    int index = Integer.parseInt(indexStr);
//    if (index != -1)
//      result.put(index, name);
//  }

  @Override
  public String rTemplateResourceURL()
  {
    return "/rejfree/plotIntensities.r";
  }

}
