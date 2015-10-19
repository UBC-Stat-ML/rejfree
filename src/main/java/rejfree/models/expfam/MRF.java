package rejfree.models.expfam;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.jblas.DoubleMatrix;
import org.jgrapht.UndirectedGraph;
import org.mvel2.templates.TemplateRuntime;

import com.google.common.base.Joiner;

import rejfree.models.normal.NormalFactor;
import bayonet.distributions.Poisson;
import bayonet.graphs.GraphUtils;
import bayonet.math.CoordinatePacker;
import bayonet.math.JBlasUtils;
import blang.annotations.DefineFactor;
import blang.variables.IntegerVariable;
import blang.variables.RealVariable;
import briefj.BriefIO;
import briefj.collections.UnorderedPair;
import briefj.opt.OptionSet;



public class MRF
{
  @OptionSet(name = "mrf")
  public final MRFOptions options;
  
  private CoordinatePacker _packer;
  public CoordinatePacker getPacker()
  {
    if (_packer == null)
      _packer = new CoordinatePacker(new int[]{options.nRows,options.nCols});
    return _packer;
  }
  
  public int nLatentVariables()
  {
    return getPacker().max;
  }
  
  public MRF()
  {
    this.options = new MRFOptions();
  }
  
  public MRF(MRFOptions options)
  {
    this.options = options;
  }
  
  /**
   * 
   * @return A modelSpec initialized at 0's
   */
  public ModelSpec newModelSpec()
  {
    return new ModelSpec();
  }
  
  public ModelSpec newModelSpecFromGenerativeProcess(Random random)
  {
    ModelSpec result = new ModelSpec();
    result.sampleFromGenerativeProcess(random);
    return result;
  }

  private DoubleMatrix localPrecision(UnorderedPair<Integer, Integer> edge)
  {
    // NOTE: if this is changed, there will also need to be changes in 
    // the Stan template
    return new DoubleMatrix(new double[][]{{options.diag, options.offDiag},{options.offDiag,options.diag}});
  }
  
  private UndirectedGraph<Integer, UnorderedPair<Integer, Integer>> _graph = null;
  public UndirectedGraph<Integer, UnorderedPair<Integer, Integer>> graph()
  {
    if (_graph == null)
      _graph = GraphUtils.grid(getPacker());
    return _graph;
  }
  
  private MultivariateNormalDistribution _mnd = null;
  private MultivariateNormalDistribution getMVNPrior()
  {
    if (_mnd == null)
      _mnd = new MultivariateNormalDistribution(new DoubleMatrix(nLatentVariables()).data, JBlasUtils.asDoubleArray(priorCovarMatrix()));
    return _mnd;
  }
  
  private DoubleMatrix _covarMatrix = null;
  private DoubleMatrix priorCovarMatrix()
  {
    if (_covarMatrix == null)
      _covarMatrix = JBlasUtils.inversePositiveMatrix(priorPrecisionMatrix());
    return _covarMatrix;
  }
  
  private DoubleMatrix _precisionMatrix = null;
  private DoubleMatrix priorPrecisionMatrix()
  {
    if (_precisionMatrix == null)
    {
      _precisionMatrix = new DoubleMatrix(nLatentVariables(), nLatentVariables());
      for (UnorderedPair<Integer, Integer> edge : graph().edgeSet())
      {
        DoubleMatrix localPrecision = localPrecision(edge);
        final int 
          v0 = edge.getFirst(),
          v1 = edge.getSecond();
        _incrementPrecision(v0, v0, localPrecision.get(0,0));
        _incrementPrecision(v1, v1, localPrecision.get(1,1));
        _incrementPrecision(v0, v1, localPrecision.get(0,1));
        _incrementPrecision(v1, v0, localPrecision.get(1,0));
      }
    }
    return _precisionMatrix;
  }
  private void _incrementPrecision(int r, int c, double value)
  {
    _precisionMatrix.put(r, c, value + _precisionMatrix.get(r, c));
  }

  public String stanModel()
  {
    String template = BriefIO.resourceToString("/rejfree/stanMRFTemplate.txt");
    return (String) TemplateRuntime.eval(template, this);
  }
  
  public String stanLatentVariableName()
  {
    return "theta";
  }
  
  public String stanObservedVariableName()
  {
    return "observations";
  }
  
  public class ModelSpec
  {
    public final List<RealVariable> latentVariables = new ArrayList<>();
    public final List<IntegerVariable> observedVariables = new ArrayList<>();
    
    @DefineFactor
    public final List<NormalFactor> priors = new ArrayList<>();
    
    @DefineFactor(onObservations = true)
    public final List<PoissonFactor> likelihoods = new ArrayList<>();
    
    /**
     * @return R-style string with topology, observation (if there is a likelihood), etc
     */
    public String getDataString()
    {
      StringBuilder result = new StringBuilder();
      
      result.append("nVertices <- " + graph().vertexSet().size() + "\n");
      result.append("nEdges <- " + graph().edgeSet().size() + "\n");
      result.append("edgeEndPoints1 <- c("  + Joiner.on(",").join(graph().edgeSet().stream().map(p -> p.getFirst()).iterator()) + ")\n");
      result.append("edgeEndPoints2 <- c("  + Joiner.on(",").join(graph().edgeSet().stream().map(p -> p.getSecond()).iterator()) + ")\n");
      result.append("diag <- " + options.diag + "\n");
      result.append("offDiag <- " + options.offDiag + "\n");
      
      if (options.hasLikelihood())
        result.append(stanObservedVariableName() + " <- c(" + Joiner.on(",").join(observedVariables.stream().map(v -> v.getIntegerValue()).iterator()) + ")\n");
      
      return result.toString();
    }
    
    public DoubleMatrix getLatentAsDoubleMatrix()
    {
      DoubleMatrix result = new DoubleMatrix(nLatentVariables());
      for (int i = 0; i < nLatentVariables(); i++)  
        result.put(i, latentVariables.get(i).getValue());
      return result;
    }

    private void sampleFromGenerativeProcess(Random random)
    {
      // generate from MVN prior
      MultivariateNormalDistribution mnd = getMVNPrior();
      mnd.reseedRandomGenerator(31*random.nextLong());
      double[] sample = mnd.sample();
      for (int i = 0; i < nLatentVariables(); i++)
        latentVariables.get(i).setValue(sample[i]);
      
      // if there is a likelihood, sample from it
      if (options.hasLikelihood())
        for (int node = 0; node < nLatentVariables(); node++)
        {
          final double rate = Math.exp(latentVariables.get(node).getValue());
          observedVariables.get(node).setValue(Poisson.generate(random, rate));
        }
    }
    
    /**
     * Assumes the nodes are indexed 0..|graph|-1
     * @param graph
     */
    private ModelSpec()
    {
      for (int i = 0; i < nLatentVariables(); i++)
      {
        latentVariables.add(new RealVariable(0.0));
        observedVariables.add(new IntegerVariable(0));
      }
      
      // likelihood: 
      if (options.hasLikelihood())
        for (Integer node : graph().vertexSet())
          likelihoods.add(new PoissonFactor(latentVariables.get(node), observedVariables.get(node)));
      
      // prior:
      for (UnorderedPair<Integer, Integer> edge : graph().edgeSet())
      {
        List<RealVariable> connected = new ArrayList<>(2);
        connected.add(latentVariables.get(edge.getFirst()));
        connected.add(latentVariables.get(edge.getSecond()));
        priors.add(new NormalFactor(localPrecision(edge), connected));
      }
    }
  }
}
