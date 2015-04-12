package rejfree.spatial;



public class SpatialUtils
{
  /**
   * Return true if the given point is contained inside the boundary.
   * See: http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
   *
   */
  public static boolean contains(double[] limits_x, double[] limits_y, double test_x, double test_y) 
  {
    int i;
    int j;
    boolean result = false;
    for (i = 0, j = limits_x.length - 1; i < limits_x.length; j = i++) 
    {
      if ((limits_y[i] > test_y) != (limits_y[j] > test_y) &&
          (test_x < (limits_x[j] - limits_x[i]) * (test_y - limits_y[i]) / (limits_y[j]-limits_y[i]) + limits_x[i])) 
        result = !result;
    }
    return result;
  }
  
  /**
   * 
   * @param c1
   * @param c2
   * @return Distance in meters between the two coordinates
   */
  public static double distance(SpatialCoordinate c1, SpatialCoordinate c2) 
  {
    final double 
      lat1 = c1.latitude,
      lat2 = c2.latitude,
      lon1 = c1.longitude,
      lon2 = c2.longitude;
    
    double theta = lon1 - lon2;
    double dist = Math.sin(deg2rad(lat1)) * Math.sin(deg2rad(lat2)) + Math.cos(deg2rad(lat1)) * Math.cos(deg2rad(lat2)) * Math.cos(deg2rad(theta));
    dist = Math.acos(dist);
    dist = rad2deg(dist);
    dist = dist * 60 * 1.1515;
    return dist * 1.609344 * 1000.0;
  }

  public static double deg2rad(double deg) 
  {
    return (deg * Math.PI / 180.0);
  }

  public static double rad2deg(double rad) 
  {
    return (rad * 180 / Math.PI);
  }

}
