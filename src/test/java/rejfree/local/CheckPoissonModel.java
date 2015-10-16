package rejfree.local;

import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;
import org.junit.Test;

import rejfree.global.GlobalRFSampler.RFSamplerOptions;
import rejfree.models.expfam.PoissonFactor;
import rejfree.models.normal.NormalFactor;
import blang.MCMCFactory;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.mcmc.Move;
import blang.mcmc.MoveFactory;
import blang.validation.CheckStationarity;
import blang.variables.IntegerVariable;
import blang.variables.RealVariable;



public class CheckPoissonModel
{
  @Test
  public void testInvar()
  {
    for (Object modelSpec : new Object[]{new TestModel(), new TestModel2()})
    {
      
      RFSamplerOptions options = new RFSamplerOptions();
      final double maxTrajectoryLen = 100.0;
      
      ProbabilityModel model = new ProbabilityModel(modelSpec);
      
      MCMCFactory factory = new MCMCFactory();
      factory.mcmcOptions.burnIn = 0;
      factory.mcmcOptions.nMCMCSweeps = 1;
      factory.mcmcOptions.CODA = false;
      factory.mcmcOptions.progressCODA = false;
      factory.mcmcOptions.thinningPeriod = 0; 
      factory.disableNodeMoveFactory();
      factory.addMoveFactory(new MoveFactory() {
        @Override
        public List<Move> build(ProbabilityModel model)
        {
          Move m = new Move() {
            @Override
            public void execute(Random rand)
            {
              LocalRFSampler sampler = new LocalRFSampler(model, options);
              sampler.iterate(rand, Integer.MAX_VALUE, maxTrajectoryLen);
            }
            @Override
            public List<?> variablesCovered()
            {
              return model.getLatentVariables();
            }
          };
          return Collections.singletonList(m);
        }
      });
      
      CheckStationarity check = new CheckStationarity();
      check.setShowSampleSummaryStats(false);
      check.check(factory.build(model), 1000, 0.01);
    }
  }
  
  class TestModel
  {
    // variables:
    RealVariable latent = new RealVariable(1);
    IntegerVariable obs = new IntegerVariable(1);
    
    // fixed hyper-parameters:
    private DoubleMatrix precision = new DoubleMatrix(new double[]{2.0});
    
    
    // factors:
    
    @DefineFactor(onObservations = true)
    public final PoissonFactor poi = new PoissonFactor(latent, obs);
    
    @DefineFactor
    public final NormalFactor nor = NormalFactor.newUnaryFactor(precision , latent);
  }
  
  class TestModel2
  {
    // variables:
    RealVariable latent1 = new RealVariable(1);
    RealVariable latent2 = new RealVariable(1);
    IntegerVariable obs1 = new IntegerVariable(1);
    IntegerVariable obs2 = new IntegerVariable(1);
    
    // fixed hyper-parameters:
    private DoubleMatrix precision = new DoubleMatrix(new double[][]{{2.0,1.5},{1.5,1.8}});
    
    // factors:
    
    @DefineFactor(onObservations = true)
    public final PoissonFactor poi1 = new PoissonFactor(latent1, obs1);
    
    @DefineFactor(onObservations = true)
    public final PoissonFactor poi2 = new PoissonFactor(latent1, obs2);
    
    @DefineFactor
    public final NormalFactor nor = NormalFactor.newBinaryFactor(precision, latent1, latent2);
  }
}
