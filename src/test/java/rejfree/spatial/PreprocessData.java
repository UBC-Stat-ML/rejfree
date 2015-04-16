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

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.base.Joiner;

import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.Indexer;
import briefj.OutputManager;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



/**
 * NB: using convention of data: (lati,long)
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class PreprocessData implements Runnable
{
  @Option
  public int streetThreshold = 2;
  
  @Option
  public File dataInput = new File("data/accidents.csv");
  
  @Option public double discretizationLen = 300;// meters
  
  /*
   * Approx coordinates for core Plateau Mont-Royal + Latin district:
   * 
   * Parc et St Jo    45.518989, -73.594328
   * Mentana et St Jo  45.529452, -73.585079
   * Notre Dame et Wolfe 45.514899, -73.549910
   * St Sulp et la commune 45.514899, -73.549910
   */
  public double [] limitPolygon_lat = new double[]{45.518989,45.529452,45.514899,45.514899};
  public double [] limitPolygon_long = new double[]{-73.594328,-73.585079,-73.549910,-73.549910};
  
  @Option public boolean useLimits = true;
  
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
      if (!keep(Double.parseDouble(lati), Double.parseDouble(longi)) || 
          street1.isEmpty() || street2.isEmpty() || 
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
    outputDiscretized(sortedIntersections);
  }
  
  private boolean keep(double lat, double longi)
  {
    if (!useLimits)
      return true;
    return SpatialUtils.contains(limitPolygon_lat, limitPolygon_long, lat, longi);
  }

  private void outputDiscretized(Map<String, List<String>> sortedIntersections)
  {
    int unkIndex = 0;
    
    OutputManager output = new OutputManager();
    output.setOutputFolder(Results.getResultFolder());
    
    Indexer<String> streetAtCornerIndexer = new Indexer<String>();
    
    for (int streetId = 0; streetId < streetIndexer.size(); streetId++)
    {
      String street = streetIndexer.i2o(streetId);
      List<String> sortedOtherStreets = sortedIntersections.get(street);
      
      String prevStreetAtCorner = null;
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
        int count = (int) countsByIntersect.getCount(UnorderedPair.of(street, otherStreet));
        
        {
          String streetAtCorner1 = street + "@" + otherStreet; streetAtCornerIndexer.addToIndex(streetAtCorner1); int streetAtCorner1Index = streetAtCornerIndexer.o2i(streetAtCorner1);
          String streetAtCorner2 = otherStreet + "@" + street; streetAtCornerIndexer.addToIndex(streetAtCorner2); int streetAtCorner2Index = streetAtCornerIndexer.o2i(streetAtCorner2);
          
          if (streetAtCorner1Index < streetAtCorner2Index)
            output.write("streetAtCornerAccidents", 
                "streetAtCorner1", streetAtCorner1,
                "streetAtCorner2", streetAtCorner2,
                "streetAtCorner1Index", streetAtCorner1Index,
                "streetAtCorner2Index", streetAtCorner2Index,
                "accidentCount", count);
          
          output.write("previousStreetAtCorners", 
              "currentStreetAtCorner", streetAtCorner1,
              "previousStreetAtCorner", (prevStreetAtCorner == null ? "unk" : prevStreetAtCorner),
              "currentStreetAtCornerIndex", streetAtCorner1Index,
              "previousStreetAtCornerIndex", (prevStreetAtCorner == null ? -1 : streetAtCornerIndexer.o2i(prevStreetAtCorner)));
              
          prevStreetAtCorner = streetAtCorner1;
        }
        
        if (distance != -1)
        {
          while (distance > discretizationLen)
          {
            String streetAtCorner = street + "@unk" + (unkIndex++);
            streetAtCornerIndexer.addToIndex(streetAtCorner);
            int streetAtCornerIndex = streetAtCornerIndexer.o2i(streetAtCorner);
            
            output.write("streetAtCornerAccidents", 
                "streetAtCorner1", streetAtCorner,
                "streetAtCorner2", "unk",
                "streetAtCorner1Index", streetAtCornerIndex,
                "streetAtCorner2Index", -1,
                "accidentCount", 0);
            
            output.write("previousStreetAtCorners", 
                "currentStreetAtCorner", streetAtCorner,
                "previousStreetAtCorner", (prevStreetAtCorner == null ? "unk" : prevStreetAtCorner),
                "currentStreetAtCornerIndex", streetAtCornerIndex,
                "previousStreetAtCornerIndex", (prevStreetAtCorner == null ? -1 : streetAtCornerIndexer.o2i(prevStreetAtCorner)));
            
            prevStreetAtCorner = streetAtCorner;
            distance -= discretizationLen;
          }
        }
      }
    }
    output.close();
    
    File rOutputFile = Results.getFileInResultFolder(R_OUTPUT_NAME);
    PrintWriter out = BriefIO.output(rOutputFile);
    addRDumpData(Results.getFileInResultFolder(GEO_OUTPUT_NAME), out, "currentStreetAtCornerIndex", "previousStreetAtCornerIndex");
    addRDumpData(Results.getFileInResultFolder(ACCIDENTS_OUTPUT_NAME), out, "streetAtCorner1Index", "streetAtCorner2Index", "accidentCount");
    out.println("nStreetAtCorners <- " + streetAtCornerIndexer.size());
    
    out.close();
  }
  
  public static final String 
    R_OUTPUT_NAME = "output-discretized.r",
    GEO_OUTPUT_NAME = "previousStreetAtCorners.csv",
    ACCIDENTS_OUTPUT_NAME = "streetAtCornerAccidents.csv";

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

  private static void addRDumpData(File csvFile, PrintWriter output, String ... fields)
  {
    Map<String,List<String>> columns = new HashMap<>();
    for (Map<String,String> row : BriefIO.readLines(csvFile).indexCSV())
      for (String key : fields)
        BriefMaps.getOrPutList(columns, key).add((row.get(key)));
    for (String key : columns.keySet())
      output.println("" + key + " <- c(" + Joiner.on(",").join(columns.get(key)) + ")");
    output.println("n" + StringUtils.capitalize(csvFile.getName().replaceAll("[.]csv$", "")) + " <- " + BriefCollections.pick(columns.values()).size());
  }

  private double distance(double[] c1, double[] c2)
  {
    return SpatialUtils.distance(SpatialCoordinate.fromArray(c1), SpatialCoordinate.fromArray(c2));
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
