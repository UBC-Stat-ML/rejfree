package rejfree.phylo;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.internal.Maps;

import rejfree.SimpleRFSampler.RFSamplerOptions;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.Move;
import blang.processing.LogDensityProcessor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.opt.InputFile;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.PhyloHMCMove;
import conifer.moves.RealVectorMHProposal;


/**
 * Test the phylogenetic MCMC moves on a phylogenetic model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestRealData implements Runnable, Processor
{
  @InputFile
  @Option(required = true, gloss = "Location of sequence files in FASTA format.")
  public File treeFile;
  
  @InputFile
  @Option(required = true, gloss = "Location of tree file in newick format")
  public File sequencesFile;
  
  @Option(gloss = "If the Rejection Free sampler should be used.")
  public boolean useRF = true;
  
  @Option(gloss = "If the Hamiltonian Monte Carlo sampler should be used.")
  public boolean useHMC = false;
  
  @Option(gloss = "If the random walk (symmetric proposal) Metropolis sampler should be used.")
  public boolean useMH = false;
  
  @Option(gloss = "If L and epsilon should be set adaptively.")
  public boolean useAdaptiveHMC = false;
  
  @Option(gloss = "Step size of the the HMC sampler")
  public double epsilon = 0.05;
  
  @Option(gloss = "Number of steps per accept reject HMC move.")
  public int L = 100;
  
  @Option(gloss = "Number of parameter moves per auxiliary variable resampling.")
  public int nItersPerPathAuxVar = 100;
  
  @OptionSet(name = "rfoptions")
  public RFSamplerOptions rfOptions = new RFSamplerOptions();
  
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
    factory.addProcessor(this);
    factory.addProcessor(new LogDensityProcessor());
    
    if (useMH)  ;
    else        factory.excludeNodeMove(RealVectorMHProposal.class);
    
    if (useRF)  factory.addNodeMove(ExpFamParameters.class, PhyloRFMove.class);
    else        ;
    
    if (useHMC) ;
    else        factory.excludeNodeMove(PhyloHMCMove.class);
    
    int nMovesRequested = (useRF ? 1 : 0) + (useHMC ? 1 : 0) + (useMH ? 1 : 0);
    
    MCMCAlgorithm mcmc = factory.build(model, false);
    if (mcmc.sampler.moves.size() != nMovesRequested)
      throw new RuntimeException("" + mcmc.sampler.moves.size() + "!=" + nMovesRequested);
    
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
        {
          ((PhyloRFMove) move).options = this.rfOptions;
          ((PhyloRFMove) move).nItersPerPathAuxVar = this.nItersPerPathAuxVar;
        }
    
    if (useHMC)
      for (Move move : mcmc.sampler.moves)
        if (move instanceof PhyloHMCMove)
          ((PhyloHMCMove) move).nItersPerPathAuxVar = this.nItersPerPathAuxVar;
    
    System.out.println(mcmc.model);
    System.out.println(mcmc.sampler);
    
    mcmc.run();
  }

  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new TestRealData());
  }
  
  Map<String,List<Double>> data = Maps.newLinkedHashMap();

  @Override
  public void process(ProcessorContext context)
  {
    
    
  }

}
