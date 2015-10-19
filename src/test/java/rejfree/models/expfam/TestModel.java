package rejfree.models.expfam;

import org.junit.Assert;
import org.junit.Test;

import bayonet.graphs.GraphUtils;
import bayonet.math.CoordinatePacker;



public class TestModel
{
  @Test
  public void test()
  {
    Assert.assertTrue(
      GraphUtils.grid(new CoordinatePacker(new int[]{3,3})).edgeSet().size() == 12 
    );
  }
}
