package rejfree.phylo;


import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;

import rejfree.SimpleRFSampler;
import rejfree.SimpleRFSampler.SimpleRFSamplerOptions;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.NodeMove;
import blang.mcmc.SampledVariable;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.ctmc.expfam.ExpectedStatistics;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.PhyloHMCMove;



public class PhyloRFMove extends NodeMove
{
  @SampledVariable ExpFamParameters parameters;
  
  @ConnectedFactor UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood;
  @ConnectedFactor IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior;
  
  public SimpleRFSamplerOptions options = new SimpleRFSamplerOptions();

  @Override
  public void execute(Random rand)
  {
    if (prior.marginalDistributionParameters.mean.getValue() != 0.0)
      throw new RuntimeException();
    final double variance = prior.marginalDistributionParameters.variance.getValue();
    
    List<PathStatistics> pathStatistics = likelihood.evolutionaryModel.samplePosteriorPaths(rand, likelihood.observations, likelihood.tree);
    
    ExpectedStatistics<CTMCState> convertedStat = PhyloHMCMove.convert(pathStatistics, parameters, likelihood); 
    CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective objective = parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0/variance, convertedStat);
    
    double [] initialPoint = parameters.getVector();
    
    
    SimpleRFSampler sampler;
    
    if (initialized)
    {
      sampler = new SimpleRFSampler(objective, new DoubleMatrix(initialPoint), options);
    }
    else
    {
      System.out.println("Initializing RF sampler");
      sampler = SimpleRFSampler.initializeRFWithLBFGS(objective, options);
      initialized = true;
    }
    sampler.iterate(rand, 1000);
    double [] newPoint = sampler.getCurrentPosition().data;
    
    parameters.setVector(newPoint);
  }

  private boolean initialized = false;


}
