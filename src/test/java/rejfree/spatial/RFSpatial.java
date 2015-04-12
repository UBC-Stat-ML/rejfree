package rejfree.spatial;

import java.util.List;

import rejfree.NormalFactor;
import blang.annotations.DefineFactor;
import briefj.run.Mains;


/**
 * Test the phylogenetic MCMC moves on a phylogenetic model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class RFSpatial implements Runnable
{

  
  public class Model
  {
    
    
//    @DefineFactor
//    public final List<NormalFactor> geographicPrior = geographicPrior();
//    
//    @DefineFactor(onObservations = true)
//    public final List<ConvolvedPoissonFactor> likelihood = likelihood();
  }
  
  public Model model;

  @Override
  public void run()
  {
    model = new Model();
    
    
  }

  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new RFSpatial());
  }


}
