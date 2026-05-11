#ifndef COMPARE_H
#define COMPARE_H

#include <cmath>
#include <string_view>

// compare
// clang-format off
constexpr auto is_not_equal = [](auto a, auto b) noexcept { return  (a < b || b < a); };
constexpr auto is_equal     = [](auto a, auto b) noexcept { return  !(a < b || b < a); };
const auto fp_is_not_equal  = [](auto a, auto b) noexcept { return ((std::isnan(a) || std::isnan(b)) ? !(std::isnan(a) && std::isnan(b)) : is_not_equal(a, b)); };
const auto fp_is_equal      = [](auto a, auto b) noexcept { return ((std::isnan(a) || std::isnan(b)) ? (std::isnan(a) && std::isnan(b)) : is_equal(a, b)); };
// clang-format on

static inline bool
cdo_cmpstr(std::string_view lhs, std::string_view rhs)
{
  return (lhs.compare(rhs) == 0);
}

#endif /* COMPARE_H */
