package rejfree.processors;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import rejfree.local.LocalRFSampler;
import rejfree.local.TrajectoryRay;
import blang.variables.RealVariable;
import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.CSV;



public class SaveRaysProcessor implements RayProcessor
{
  private final Map<RealVariable, List<TrajectoryRay>> samples;
  private double totalLength = 0.0;
  
  public SaveRaysProcessor(Collection<RealVariable> variables)
  {
    if (variables.size() == 1) // slight efficiency speed-up
      this.samples = Collections.singletonMap(BriefCollections.pick(variables), new ArrayList<>());
    else
    {
      this.samples = new LinkedHashMap<>();
      for (RealVariable key : variables)
        samples.put(key, new ArrayList<>());
    }
  }
  
  public void toCSV(File f)
  {
    PrintWriter output = BriefIO.output(f);
    output.println(CSV.toCSV("variableIndex", "time", "position"));
    int varIndex = 0;
    for (List<TrajectoryRay> traj : samples.values())
    {
      for (TrajectoryRay ray : traj)
        output.println(CSV.toCSV(varIndex, ray.t, ray.position_t));
      varIndex++;
    }
    output.close();
  }

  @Override
  public void init(LocalRFSampler sampler)
  {
  }
  
  public double totalLength()
  {
    return totalLength;
  }
  
  public List<Double> convertToSample(RealVariable var, double delta)
  {
    return convertToSamples(samples.get(var), delta, totalLength);
  }
  
  private static List<Double> convertToSamples(List<TrajectoryRay> rays, double delta, double totalLength)
  {
    List<Double> result = new ArrayList<>();
    double timeToNextCollect = 0.0;
    for (int i = 0; i < rays.size(); i++)  
    {
      TrajectoryRay 
        curRay = rays.get(i),
        nxtRay = i == rays.size() - 1 ? null : rays.get(i+1);
      double endTimeForRay = (nxtRay == null ? totalLength : nxtRay.t);
      for (double t = curRay.t + timeToNextCollect; t < endTimeForRay; t += delta)
        result.add(curRay.position(t));
      double rayLen = endTimeForRay - curRay.t;
      if (rayLen < timeToNextCollect)
        timeToNextCollect -= rayLen;
      else
        timeToNextCollect = rayLen - delta * ((int) (rayLen / delta));
    }
    return result;
  }
  
  @Override
  public void processRay(RealVariable var, TrajectoryRay ray, double time,
      LocalRFSampler sampler)
  {
    List<TrajectoryRay> list = samples.get(var);
    if (list != null)
      list.add(ray);
    this.totalLength = time;
  }
  
}