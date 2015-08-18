package rejfree.models.normal;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.jblas.DoubleMatrix;
import org.mvel2.templates.TemplateRuntime;

import rejfree.local.CollisionFactor;
import bayonet.bugs.StanWrapper;
import bayonet.math.JBlasUtils;
import blang.annotations.DefineFactor;
import blang.variables.RealVariable;
import briefj.BriefIO;



public class NormalChain
{
  public final NormalChainOptions options;
  
  public NormalChain(NormalChainOptions options)
  {
    this.options = options;
    buildPrecisionMatrices();
  }

  public int dim() { return options.nPairs + 1; }
  
  private DoubleMatrix fullPrecision;
  public DoubleMatrix covarMatrix;
  private List<DoubleMatrix> pairPrecisions;
  private MultivariateNormalDistribution normal;
  
  public class NormalChainModel
  {
    public List<RealVariable> variables = new ArrayList<>();
    
    @DefineFactor
    public final List<CollisionFactor> localFactors;
    
    public NormalChainModel(double [] init)
    {
      this.localFactors = options.useLocal ?
        localFactors(init) :
        globalFactor(init);  
    }
    
    private List<CollisionFactor> globalFactor(double [] init)
    {
      List<CollisionFactor> result = new ArrayList<>();
      
      for (int i = 0; i < init.length; i++)
      {
        RealVariable current = RealVariable.real(init[i]);
        variables.add(current);
      }
      
      CollisionFactor f = 
          options.useAnalytic ? 
          new NormalFactor(fullPrecision, variables) :   
          NumericNormalFactor.withPrecision(variables, fullPrecision);
      
      result.add(f);
      
      return result;
    }

    private List<CollisionFactor> localFactors(double [] init)
    {
      List<CollisionFactor> result = new ArrayList<>();
      RealVariable prev = RealVariable.real(init[0]);
      variables.add(prev);
      for (int i = 0; i < init.length - 1; i++)
      {
        DoubleMatrix localPrec = pairPrecisions.get(i);
        
        RealVariable current = RealVariable.real(init[i+1]);
        variables.add(current);
        
        List<RealVariable> variables = new ArrayList<>();
        variables.add(prev);
        variables.add(current);
        
        CollisionFactor f = 
            options.useAnalytic ? 
            NormalFactor.newBinaryFactor(localPrec, current, prev) :
            NumericNormalFactor.withPrecision(variables, localPrec);
        
        result.add(f);
        prev = current;
      }
      
      return result;
    }
  }
  
  public DoubleMatrix exactSample()
  {
    return new DoubleMatrix(normal.sample());
  }
  
  private void buildPrecisionMatrices()
  {
    fullPrecision = new DoubleMatrix(options.nPairs+1, options.nPairs+1);
    pairPrecisions = new ArrayList<DoubleMatrix>();
    for (int i = 0; i < options.nPairs; i++)
    {
      DoubleMatrix cur = new DoubleMatrix(2,2);
      pairPrecisions.add(cur);
      for (int r = 0; r < 2; r++)
        for (int c = 0; c < 2; c++)
        {
          double val = r == c ? options.diag : options.offDiag;
          cur.put(r, c, val);
          fullPrecision.put(r+i,c+i,val + fullPrecision.get(r+i,c+i));
        }
    }
    covarMatrix = JBlasUtils.inversePositiveMatrix(fullPrecision);  
    normal = new MultivariateNormalDistribution(new DoubleMatrix(dim()).data, JBlasUtils.asDoubleArray(covarMatrix));
    normal.reseedRandomGenerator(options.random.nextLong());
  }
  
  private String stanModel()
  {
    String template = BriefIO.resourceToString("/rejfree/stanChainTemplate.txt");
    return (String) TemplateRuntime.eval(template, this);
  }
  
  public File stanProgram(File stanHome)
  {
    return StanWrapper.compile(stanModel(), stanHome);
  }

}