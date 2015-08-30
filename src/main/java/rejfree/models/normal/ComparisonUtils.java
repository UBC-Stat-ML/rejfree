package rejfree.models.normal;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;



public class ComparisonUtils
{
  public static <T> Collection<T> subsample(Collection<T> items, int interval)
  {
    List<T> result = new ArrayList<>();
    int current = 0;
    for (T item : items)
      if (current++ % interval == 0)
        result.add(item);
    return result;
  }
}
