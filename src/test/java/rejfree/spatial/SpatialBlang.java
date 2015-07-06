package rejfree.spatial;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.local.LocalRFSampler;
import rejfree.local.NormalFactor;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.IntegerVariable;
import blang.variables.RealVariable;
import briefj.BriefIO;
import briefj.opt.OptionSet;
import briefj.run.Results;

import com.google.common.base.Joiner;


/**
 * Test the phylogenetic MCMC moves on a phylogenetic model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class SpatialBlang implements Processor
{
  public final SpatialMainOptions mainOptions;
  
  public SpatialBlang(SpatialMainOptions mainOptions)
  {
    this.mainOptions = mainOptions;
  }

  @OptionSet(name = "rf")
  public final RFSamplerOptions options = new RFSamplerOptions();
  
  public class Model
  {
    // indexed by streetAtCornerIndex in the csv files
    Map<Integer, RealVariable> variables = variables();
    
    @DefineFactor
    public final List<NormalFactor> geographicPrior = geographicPrior(mainOptions.getPriorTransitionVariance(), mainOptions.getPriorInitialVariance());
    
    @DefineFactor(onObservations = true)
    public final List<ConvolvedPoissonFactor> likelihood = likelihood();

    private Map<Integer, RealVariable> variables()
    {
      Map<Integer, RealVariable> result = new LinkedHashMap<Integer, RealVariable>();
      for (Map<String,String> line : BriefIO.readLines(mainOptions.getGeoDataCSVFile()).indexCSV())
        result.put(Integer.parseInt(line.get("currentStreetAtCornerIndex")), new RealVariable(0.0));
      return result;
    }

    private List<ConvolvedPoissonFactor> likelihood()
    {
      List<ConvolvedPoissonFactor> result = new ArrayList<>();
      
      for (Map<String,String> line : BriefIO.readLines(mainOptions.getAccidentsDataCSVFile()).indexCSV())
      {
        RealVariable 
          first  = variables.get(Integer.parseInt(line.get("streetAtCorner1Index"))),
          second = variables.get(Integer.parseInt(line.get("streetAtCorner2Index")));
        
        int datum = Integer.parseInt(line.get("accidentCount"));
        
        ConvolvedPoissonFactor factor = new ConvolvedPoissonFactor(first, second, new IntegerVariable(datum));
        result.add(factor);
      }
      
      return result;
    }

    @SuppressWarnings("unused")
    private List<NormalFactor> geographicPrior(double priorTransitionVariance, double priorInitialVariance)
    {
      List<NormalFactor> result = new ArrayList<>();
      
      for (Map<String,String> line : BriefIO.readLines(mainOptions.getGeoDataCSVFile()).indexCSV())
      {
        RealVariable 
          current = variables.get(Integer.parseInt(line.get("currentStreetAtCornerIndex"))),
          prev    = variables.get(Integer.parseInt(line.get("previousStreetAtCornerIndex")));
        
//        List<RealVariable> variables;
//        double [][] covar;
//        if (prev == null)
//        {
//          // beginning of a chain
//          variables = Collections.singletonList(current);
//          covar = new double[][]{{init}};
//        }
//        else
//        {
//          // continuation
//          variables = new ArrayList<>();
//          variables.add(current);
//          variables.add(prev);
//          // TODO: fix this hack!
//          double entry = 1.0/variance;
//          covar = new double[][]{{entry,-entry},{-entry,entry}};
//        }
//        NormalFactor f = NormalFactor.withPrecision(variables, new DoubleMatrix(covar));
        
        NormalFactor f = null; if (true) throw new RuntimeException("Should fix and uncomment code below");
//        prev == null ? 
//          SpatialNormalFactor.newUnaryFactor(priorInitialVariance, current) : 
//          SpatialNormalFactor.newBinaryFactor(priorTransitionVariance, current, prev);
        
        result.add(f);
      }
      
      return result;
    }
  }
  
  public Model modelSpec;
  private PrintWriter output;

  public void runRF()
  {
    ProbabilityModel model = setupModel();
    
    LocalRFSampler local = new LocalRFSampler(model, options);
    local.addPointProcessor(this);
    local.iterate(mainOptions.random, mainOptions.nSamples);
    
    tearDown();
  }
  
  public void runStandardBlang()
  {
    ProbabilityModel model = setupModel();
    
    MCMCFactory factory = mainOptions.getMCMCFactory();
    factory.addProcessor(this);
    MCMCAlgorithm mcmcAlgo = factory.build(model);
    mcmcAlgo.run();
    
    tearDown();
  }
  
  private void tearDown()
  {
    output.close();
  }
  
  private ProbabilityModel setupModel()
  {
    modelSpec = new Model();
    ProbabilityModel model = new ProbabilityModel(modelSpec);
    System.out.println(model);
    File samplesFile = Results.getFileInResultFolder(RunSpatialExample.SAMPLES_FILE_NAME);
    output = BriefIO.output(samplesFile);
    printHeader(output);
    return model;
  }

  private void printHeader(PrintWriter output)
  {
    output.println(Joiner.on(",").join(
        IntStream.range(0, modelSpec.variables.size()).
          mapToObj(i -> "logIntensity." + (i+1)).
          collect(Collectors.toList())));
  }

  @Override
  public void process(ProcessorContext context)
  {
    output.println(Joiner.on(",").join(
        IntStream.range(0, modelSpec.variables.size()).
          mapToObj(i -> modelSpec.variables.get(i).getValue()).
          collect(Collectors.toList())));
  }
}
