#ifndef CDO_MATH_H
#define CDO_MATH_H

#include <array>

// clang-format off
namespace cdo
{
  constexpr double sqr(double x) noexcept { return x * x; }
  constexpr double sqr_distance(std::array<double, 3> const &a, std::array<double, 3> const &b) noexcept
  {
    return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
  }
  constexpr double sqr_distance(double const (&a)[3], double const (&b)[3]) noexcept
  {
    return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
  }
  unsigned int is_power_of_two(unsigned int  x);
  double NaN();
}
// clang-format on

#endif
