#ifndef GRID_CONVERT_H
#define GRID_CONVERT_H

#include <cmath>
#include <numbers>
#include <array>

static inline void
gcLLtoXYZ(double lon, double lat, std::array<double, 3> &xyz)
{
  auto cos_lat = std::cos(lat);
  xyz[0] = cos_lat * std::cos(lon);
  xyz[1] = cos_lat * std::sin(lon);
  xyz[2] = std::sin(lat);
}

static inline void
gcLLtoXYZ(double lon, double lat, double (&xyz)[3])
{
  auto cos_lat = std::cos(lat);
  xyz[0] = cos_lat * std::cos(lon);
  xyz[1] = cos_lat * std::sin(lon);
  xyz[2] = std::sin(lat);
}

#ifndef RAD_CONVERT
#define RAD_CONVERT
constexpr double RAD2DEG = 180.0 / std::numbers::pi;  // conversion from radians to degrees
constexpr double DEG2RAD = std::numbers::pi / 180.0;  // conversion from degrees to radians
// clang-format off
inline constexpr double rad_to_deg(double rad) noexcept { return rad * 180.0 / std::numbers::pi; };
inline constexpr double deg_to_rad(double deg) noexcept { return deg * std::numbers::pi / 180.0; };
// clang-format on
#endif

#endif
