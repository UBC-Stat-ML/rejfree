package rejfree.phylo;

import java.io.File;

import rejfree.SimpleRFSampler.SimpleRFSamplerOptions;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.MHMove;
import blang.mcmc.Move;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.PhyloHMCMove;


/**
 * Test the phylogenetic MCMC moves on a simple tree model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestRealData implements Runnable
{
  @Option 
  public File treeFile = new File("primates.fasta");
  
  @Option
  public File sequencesFile = new File("final-tree.newick");
  
  @Option 
  public boolean useRF = false;
  
  @Option 
  public boolean useHMC = false;
  
  @Option
  public boolean useAdaptiveHMC = false;
  
  @Option
  public double epsilon = 0.1;
  
  @Option
  public int L = 10;
  
  @OptionSet(name = "rfoptions")
  SimpleRFSamplerOptions rfOptions = new SimpleRFSamplerOptions();
  
  @OptionSet(name = "factory")
  public final MCMCFactory factory = new MCMCFactory();
  
  public class Model
  {
    @DefineFactor(onObservations = true)
    public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
      UnrootedTreeLikelihood.fromFastaFile(treeFile).withExpFamMixture(ExpFamMixture.dnaGTR()).withTree(sequencesFile);
    
    @DefineFactor
    public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> priorOnParams =
      IIDRealVectorGenerativeFactor.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
  }
  
  public Model model;

  @Override
  public void run()
  {
    model = new Model();
    factory.excludeNodeMove(MHMove.class);
    if (useRF)
      factory.addNodeMove(ExpFamParameters.class, PhyloRFMove.class);
    if (useHMC)
      ;
    else
      factory.excludeNodeMove(PhyloHMCMove.class);
    int nMovesRequested = (useRF ? 1 : 0) + (useHMC ? 1 : 0);
    
    MCMCAlgorithm mcmc = factory.build(model, false);
    if (mcmc.sampler.moves.size() != nMovesRequested)
      throw new RuntimeException();
    
    if (!useAdaptiveHMC && useHMC)
    {
      // find the hmc sampler
      PhyloHMCMove hmcMove = null;
      for (Move move : mcmc.sampler.moves)
        if (move instanceof PhyloHMCMove)
          hmcMove = (PhyloHMCMove) move;
      // disable adaptivity by setting fixed L, epsilon values
      hmcMove.epsilon = this.epsilon;
      hmcMove.L = this.L;
    }
    
    if (useRF)
      for (Move move : mcmc.sampler.moves)
        if (move instanceof PhyloRFMove)
          ((PhyloRFMove) move).options = this.rfOptions;
    
    System.out.println(mcmc.model);
    System.out.println(mcmc.sampler);
    mcmc.run();
  }

  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new TestRealData());
  }

}
