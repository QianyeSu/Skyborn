#ifndef POINT_H
#define POINT_H

class PointLonLat
{
public:
  PointLonLat() {};
  PointLonLat(double lon, double lat) : m_lon{ lon }, m_lat{ lat } {};

  // clang-format off
  double lon() const noexcept { return m_lon; };
  double lat() const noexcept { return m_lat; };
  // clang-format on

private:
  double m_lon{ 0.0 };
  double m_lat{ 0.0 };
};

#endif
