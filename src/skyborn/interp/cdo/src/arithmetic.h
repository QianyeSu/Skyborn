#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include <cmath>

// clang-format off

namespace arithmetic
{
const auto add = [](double x, double y) noexcept { return x + y; };
const auto sub = [](double x, double y) noexcept { return x - y; };
const auto mul = [](double x, double y) noexcept { return x * y; };
const auto div = [](double x, double y) noexcept { return x / y; };
const auto pow = [](double x, double y) noexcept { return std::pow(x, y); };
const auto sqrt = [](double x) noexcept { return std::sqrt(x); };

const auto addm = [](auto x, auto y, auto mvx, auto mvy, auto is_EQ) { return (is_EQ(x, mvx) || is_EQ(y, mvy)) ? mvx : x + y; };
const auto subm = [](auto x, auto y, auto mvx, auto mvy, auto is_EQ) { return (is_EQ(x, mvx) || is_EQ(y, mvy)) ? mvx : x - y; };
const auto mulm = [](auto x, auto y, auto mvx, auto mvy, auto is_EQ) { return (is_EQ(x, 0) || is_EQ(y, 0)) ? 0 : (is_EQ(x, mvx) || is_EQ(y, mvy)) ? mvx : x * y; };
const auto divm = [](auto x, auto y, auto mvx, auto mvy, auto is_EQ) { return (is_EQ(x, mvx) || is_EQ(y, mvy) || is_EQ(y, 0)) ? mvx : x / y; };
const auto powm = [](auto x, auto y, auto mvx, auto mvy, auto is_EQ) { return (is_EQ(x, mvx) || is_EQ(y, mvy)) ? mvx : std::pow(x, y); };
const auto sqrtm = [](auto x, auto mvx, auto is_EQ) { return (is_EQ(x, mvx) || x < 0) ? mvx : std::sqrt(x); };
const auto divmx = [](auto x, auto y, auto mvx, auto is_EQ) { return (is_EQ(x, mvx) || is_EQ(y, 0)) ? mvx : x / y; };
const auto divx = [](auto x, auto y, auto mvx, auto is_EQ) { return is_EQ(y, 0) ? mvx : x / y; };
}

const auto ADD = arithmetic::add;
const auto SUB = arithmetic::sub;
const auto MUL = arithmetic::mul;
const auto DIV = arithmetic::div;
const auto POW = arithmetic::pow;
const auto SQRT = arithmetic::sqrt;

#define ADDM(x, y)  arithmetic::addm(x, y, missval1, missval2, is_EQ)
#define SUBM(x, y)  arithmetic::subm(x, y, missval1, missval2, is_EQ)
#define MULM(x, y)  arithmetic::mulm(x, y, missval1, missval2, is_EQ)
#define DIVM(x, y)  arithmetic::divm(x, y, missval1, missval2, is_EQ)
#define POWM(x, y)  arithmetic::powm(x, y, missval1, missval2, is_EQ)
#define SQRTM(x)    arithmetic::sqrtm(x, missval1, is_EQ)
#define DIVMX(x, y)  arithmetic::divmx(x, y, missval1, is_EQ)
#define DIVX(x, y)  arithmetic::divx(x, y, missval1, is_EQ)

// clang-format on

#endif
