package rejfree;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;



public class EventQueue<S>
{
  private final TreeMap<Double,S> sortedEvents = new TreeMap<>();
  private final Map<S, Double> eventTimes = new HashMap<>();
  
  public Entry<Double,S> pollEvent()
  {
    return sortedEvents.pollFirstEntry();
  }
  
  public void remove(S event)
  {
    double time = eventTimes.get(event);
    sortedEvents.remove(time);
    eventTimes.remove(event);
  }
  
  public void add(S event, double time)
  {
    if (sortedEvents.containsKey(time))
      throw new RuntimeException("EventQueue does not support two events at the same time (t=" + time + ",event=" + event + ")");
    sortedEvents.put(time, event);
    eventTimes.put(event, time);
  }

  public double peekTime()
  {
    return sortedEvents.firstKey();
  }
}