package rejfree.processors;

import java.util.ArrayList;
import java.util.List;

import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.RealVariable;



public class SaveSamplesProcessor implements Processor
{
  private final List<RealVariable> variables;
  public final List<List<Double>> samples;
  private final int nVars;
  
  public SaveSamplesProcessor(List<RealVariable> variables)
  {
    this.variables = variables;
    this.nVars = variables.size();
    this.samples = new ArrayList<>();
    for (int varIndex = 0; varIndex < nVars; varIndex++)
      this.samples.add(new ArrayList<>());
  }

  @Override
  public void process(ProcessorContext context)
  {
    for (int varIndex = 0; varIndex < nVars; varIndex++)
      samples.get(varIndex).add(variables.get(varIndex).getValue());
  }
  
}