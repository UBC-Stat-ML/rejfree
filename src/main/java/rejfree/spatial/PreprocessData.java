package rejfree.spatial;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;

import com.google.common.base.Joiner;

import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class PreprocessData implements Runnable
{
  @Option
  public int streetThreshold = 10;
  
  @Option
  public File dataInput = new File("data/accidents.csv");
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new PreprocessData());
  }
  
  private final Counter<UnorderedPair<String, String>> countsByIntersect = new Counter<>();
  private final Map<UnorderedPair<String, String>, Pair<Double,Double>> cornerLocations = new LinkedHashMap<>();
  private final Map<String,Set<String>> intersectionsMap = new HashMap<>();
  private final Indexer<String> streetIndexer = new Indexer<>();
  
  private void loadData()
  {
    Counter<String> countsByStreet = getCountByStreet();
    loop : for (Map<String,String> fields : BriefIO.readLines(dataInput).indexCSV())
    {
      String 
        street1 = fields.get("rue1"),
        street2 = fields.get("rue2"),
        lati  = fields.get("lat"),
        longi = fields.get("long");
      street1 = cleanStreetName(street1);
      street2 = cleanStreetName(street2);
      if (street1.isEmpty() || street2.isEmpty() || 
          countsByStreet.getCount(street1) < streetThreshold ||
          countsByStreet.getCount(street2) < streetThreshold)
        continue loop;
      
      UnorderedPair<String, String> key = new UnorderedPair<String, String>(street1, street2);
      countsByIntersect.incrementCount(key, 1.0);
      Pair<Double,Double> location = Pair.of(Double.parseDouble(lati), Double.parseDouble(longi));
      // warn for contradictory info in dataset
      if (cornerLocations.containsKey(key) && !cornerLocations.get(key).equals(location))
        System.err.println("Conflicting spatial info for corner of " + key + ":" 
              + cornerLocations.get(key) + " vs " + location
              + "\n\tKeeping the latter.");
      cornerLocations.put(key, location);
      
      BriefMaps.getOrPutSet(intersectionsMap, street1).add(street2);
      BriefMaps.getOrPutSet(intersectionsMap, street2).add(street1);
      
      streetIndexer.addToIndex(street1);
      streetIndexer.addToIndex(street2);
    }
    System.out.println("counstByStreet: " + countsByStreet);
  }
  
  public void run()
  {
    loadData();
    Map<String,List<String>> sortedIntersections = sortIntersections();
    output(sortedIntersections);
  }

  private void output(Map<String, List<String>> sortedIntersections)
  {
    File csvOutputFile = Results.getFileInResultFolder("output.csv");
    PrintWriter 
      output  = BriefIO.output(csvOutputFile);
    
    output.println("streetId,streetName,intersectionIndexForThisStreet,otherStreetName,otherStreetIndex,distanceToNextInThisStreet,count");
    for (int streetId = 0; streetId < streetIndexer.size(); streetId++)
    {
      String street = streetIndexer.i2o(streetId);
      List<String> sortedOtherStreets = sortedIntersections.get(street);
      
      for (int intersectionIndexForThisStreet = 0; intersectionIndexForThisStreet < sortedOtherStreets.size(); intersectionIndexForThisStreet++)
      {
        String otherStreet = sortedOtherStreets.get(intersectionIndexForThisStreet);
        double distance = -1.0;
        if (intersectionIndexForThisStreet < sortedOtherStreets.size() - 1)
        {
          String nextStreet = sortedOtherStreets.get(intersectionIndexForThisStreet + 1);
          distance = distance(
              asArray(cornerLocations.get(UnorderedPair.of(street, otherStreet))),
              asArray(cornerLocations.get(UnorderedPair.of(street, nextStreet))));
        }
        int otherStreetId = streetIndexer.o2i(otherStreet);
        int count = (int) countsByIntersect.getCount(UnorderedPair.of(street, otherStreet));
        output.println(Joiner.on(",").join(streetId, street, intersectionIndexForThisStreet, otherStreet, otherStreetId, distance, count));
      }
    }
    output.close();
    File rOutputFile = Results.getFileInResultFolder("output.r");
    createRDumpFile(csvOutputFile, rOutputFile);
  }

  private Map<String, List<String>> sortIntersections()
  {
    Map<String, List<String>> sortedIntersections = new LinkedHashMap<>();
    for (String street : streetIndexer.objects())
    {
      // order the intersections within each street
      PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
      pca.setup(intersectionsMap.get(street).size(), 2);
      for (String otherStreet : intersectionsMap.get(street))
      {
        Pair<Double,Double> coordinate = cornerLocations.get(UnorderedPair.of(street, otherStreet));
        double [] vector = new double[]{coordinate.getLeft(), coordinate.getRight()};
        pca.addSample(vector);
      }
      pca.computeBasis(1);
      
      List<Pair<String,Double>> sorted = new ArrayList<>();
      for (String otherStreet : intersectionsMap.get(street))
      {
        Pair<Double,Double> coordinate = cornerLocations.get(UnorderedPair.of(street, otherStreet));
        double [] vector = new double[]{coordinate.getLeft(), coordinate.getRight()};
        sorted.add(Pair.of(otherStreet, pca.sampleToEigenSpace(vector)[0]));
      }
      Collections.sort(sorted, (item1, item2) -> item1.getRight().compareTo(item2.getRight()));
      
      System.out.println(street + " : " );
      Pair<Double,Double> prev = null;
      List<String> sortedNames = new ArrayList<>();
      for (Pair<String,Double> other : sorted)
      {
        Pair<Double,Double> cur = cornerLocations.get(UnorderedPair.of(other.getLeft(), street));
        System.out.println("\t" + other.getLeft() + "\t" + delta(prev, cur));
        prev = cur;
        sortedNames.add(other.getLeft());
      }
      sortedIntersections.put(street, sortedNames);
    }
    return sortedIntersections;
  }

  private Counter<String> getCountByStreet()
  {
    Counter<String> countsByStreet = new Counter<>();
    for (Map<String,String> fields : BriefIO.readLines(dataInput ).indexCSV())
    {
      String 
        street1 = fields.get("rue1"),
        street2 = fields.get("rue2");
      street1 = cleanStreetName(street1);
      street2 = cleanStreetName(street2);
      countsByStreet.incrementCount(street1, 1.0);
      countsByStreet.incrementCount(street2, 1.0);
    }
    return countsByStreet;
  }

  private static void createRDumpFile(File csvOutputFile, File rOutputFile)
  {
    Map<String,List<String>> columns = new HashMap<>();
    for (Map<String,String> row : BriefIO.readLines(csvOutputFile).indexCSV())
      for (String key : new String [] {"otherStreetIndex", "intersectionIndexForThisStreet", "count", "distanceToNextInThisStreet", "streetId"})
        BriefMaps.getOrPutList(columns, key).add((row.get(key)));
    PrintWriter output = BriefIO.output(rOutputFile);
    for (String key : columns.keySet())
      output.println("" + key + " <- c(" + Joiner.on(",").join(columns.get(key)) + ")");
    output.println("N <- " + BriefCollections.pick(columns.values()).size());
    output.close();
  }

  private double distance(double[] asArray, double[] asArray2)
  {
    DoubleMatrix v1 = new DoubleMatrix(asArray), v2 = new DoubleMatrix(asArray2);
    return v1.distance2(v2);
  }

  private static String delta(
      Pair<Double, Double> prev,
      Pair<Double, Double> cur)
  {
    if (prev == null)
      return "";
    return sign(prev, cur, 0) + " " + sign(prev, cur, 1);
  }

  private static String sign(Pair<Double, Double> prev, Pair<Double, Double> cur,
      int i)
  {
    return asArray(prev)[i] < asArray(cur)[i] ? "+" : "-";
  }
  
  private static double [] asArray(Pair<Double,Double> item) 
  {
    return new double[]{item.getLeft(), item.getRight()};
  }

  private static String cleanStreetName(String string)
  {
    if (string == null)
      return null;
    string = string.replaceAll("\\s*$", "");
    string = string.replaceAll("^\\s*", "");
    return string.replaceAll(" [NSEO]$", "").toUpperCase().replaceAll("\\s+", "_");
  }
}
