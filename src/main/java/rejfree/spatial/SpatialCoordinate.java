package rejfree.spatial;



public class SpatialCoordinate
{
  public final double latitude, longitude;

  private SpatialCoordinate(double latitude, double longitude)
  {
    if (latitude < -90.0 || latitude > 90 || longitude < -180 || longitude > 180 )
      throw new RuntimeException();
    this.latitude = latitude;
    this.longitude = longitude;
  }
  
  public static SpatialCoordinate fromArray(double [] array)
  {
    if (array.length != 2)
      throw new RuntimeException();
    return new SpatialCoordinate(array[0], array[1]);
  }

  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    long temp;
    temp = Double.doubleToLongBits(latitude);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    temp = Double.doubleToLongBits(longitude);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    return result;
  }

  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    SpatialCoordinate other = (SpatialCoordinate) obj;
    if (Double.doubleToLongBits(latitude) != Double
        .doubleToLongBits(other.latitude))
      return false;
    if (Double.doubleToLongBits(longitude) != Double
        .doubleToLongBits(other.longitude))
      return false;
    return true;
  }

  @Override
  public String toString()
  {
    return "SpatialCoordinate [latitude=" + latitude + ", longitude="
        + longitude + "]";
  }
}