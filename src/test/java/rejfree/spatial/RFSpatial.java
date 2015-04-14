package rejfree.spatial;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.jblas.DoubleMatrix;

import rejfree.GlobalRFSampler.RFSamplerOptions;
import rejfree.NormalFactor;
import rejfree.local.LocalRFSampler;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.variables.IntegerVariable;
import blang.variables.RealVariable;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;


/**
 * Test the phylogenetic MCMC moves on a phylogenetic model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class RFSpatial implements Runnable
{
  
  @Option(required = true)
  public File accidentsData;

  @Option(required = true)
  public File geographicData;
  
  @Option
  public double drift = 0.2;
  
  @Option
  public double init = 10;
  
  @Option
  public int nIterations = 1000;
  
  @Option
  public Random random = new Random(1);
  
  @OptionSet(name = "rf")
  public RFSamplerOptions options = new RFSamplerOptions();

  public class Model
  {
    // indexed by streetAtCornerIndex in the csv files
    Map<Integer, RealVariable> variables = variables();
    
    @DefineFactor
    public final List<NormalFactor> geographicPrior = geographicPrior(drift, init);
    
    @DefineFactor(onObservations = true)
    public final List<ConvolvedPoissonFactor> likelihood = likelihood();

    private Map<Integer, RealVariable> variables()
    {
      Map<Integer, RealVariable> result = new LinkedHashMap<Integer, RealVariable>();
      for (Map<String,String> line : BriefIO.readLines(geographicData).indexCSV())
        result.put(Integer.parseInt(line.get("currentStreetAtCornerIndex")), new RealVariable(0.0));
      return result;
    }

    private List<ConvolvedPoissonFactor> likelihood()
    {
      List<ConvolvedPoissonFactor> result = new ArrayList<>();
      
      for (Map<String,String> line : BriefIO.readLines(accidentsData).indexCSV())
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

    private List<NormalFactor> geographicPrior(double drift, double init)
    {
      List<NormalFactor> result = new ArrayList<NormalFactor>();
      
      for (Map<String,String> line : BriefIO.readLines(geographicData).indexCSV())
      {
        RealVariable 
          current = variables.get(Integer.parseInt(line.get("currentStreetAtCornerIndex"))),
          prev    = variables.get(Integer.parseInt(line.get("previousStreetAtCornerIndex")));
        
        List<RealVariable> variables;
        double [][] covar;
        if (prev == null)
        {
          // beginning of a chain
          variables = Collections.singletonList(current);
          covar = new double[][]{{init}};
        }
        else
        {
          // continuation
          variables = new ArrayList<>();
          variables.add(current);
          variables.add(prev);
          double entry = 1.0/drift;
          covar = new double[][]{{entry,-entry},{-entry,entry}};
        }
        NormalFactor f = NormalFactor.withPrecision(variables, new DoubleMatrix(covar));
        result.add(f);
      }
      
      return result;
    }
  }
  
  public Model modelSpec;

  @Override
  public void run()
  {
    modelSpec = new Model();
    ProbabilityModel model = new ProbabilityModel(modelSpec);
//    System.out.println(model);
    LocalRFSampler local = new LocalRFSampler(model, options);
    local.iterate(random, nIterations);
  }

  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new RFSpatial());
//    double point0 = 0.34;
//    double point1 = 0.12;
//    double drift = 1.0;
////    double init = 2.0;
//    
//    // method 1
////    double method1 = Normal.logDensity(point1, point0, drift);
//    double method1 = - 0.5 * (point0 - point1) * (point0 - point1) / drift;
//    System.out.println(method1);
//    
//    double id = 1.0/drift;
//    DoubleMatrix precision = new DoubleMatrix(new double[][]{{id, -id},{-id, id}});
//    NormalEnergy e = NormalEnergy.withPrecision(precision);
//    System.out.println(Arrays.toString(e.derivativeAt(new double[]{point0, point1})));
//    double method2 = - e.valueAt(new double[]{point0, point1});
//    
//    System.out.println(method2);
//    
  }


}
