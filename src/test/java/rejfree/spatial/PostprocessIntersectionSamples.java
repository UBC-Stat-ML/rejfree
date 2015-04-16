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
  
  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new PostprocessIntersectionSamples());
  }
  
  @Override
  public void run()
  {
    RUtils.callRBridge(this);
  }

  @Override
  public String rTemplateResourceURL()
  {
    return "/rejfree/plotIntensities.r";
  }
}
