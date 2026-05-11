/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef VARRAY_H
#define VARRAY_H

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <iostream>  // cerr
#include <limits>
#include <span>
#include <vector>
#include "cpp_lib.h"
#include "compare.h"

// unary operators
// clang-format off
const auto unary_op_abs   = [](double x) noexcept { return std::fabs(x); };
const auto unary_op_int   = [](auto   x) noexcept { return (int)(x); };
const auto unary_op_nint  = [](double x) noexcept { return std::round(x); };
const auto unary_op_sqr   = [](double x) noexcept { return x * x; };
const auto unary_op_nop   = [](auto   x) noexcept { return x; };
const auto unary_op_reci  = [](double x) noexcept { return 1.0 / x; };
const auto unary_op_not   = [](auto   x) noexcept { return is_equal(x, 0); };
const auto unary_op_exp   = [](double x) noexcept { return std::exp(x); };
const auto unary_op_log   = [](double x) noexcept { return std::log(x); };
const auto unary_op_log10 = [](double x) noexcept { return std::log10(x); };
const auto unary_op_sin   = [](double x) noexcept { return std::sin(x); };
const auto unary_op_cos   = [](double x) noexcept { return std::cos(x); };
const auto unary_op_tan   = [](double x) noexcept { return std::tan(x); };
const auto unary_op_asin  = [](double x) noexcept { return std::asin(x); };
const auto unary_op_acos  = [](double x) noexcept { return std::acos(x); };
const auto unary_op_atan  = [](double x) noexcept { return std::atan(x); };
// clang-format on

// binary operators
// clang-format off
const auto binary_op_LT  = [](auto x, auto y) noexcept { return x < y; };
const auto binary_op_GT  = [](auto x, auto y) noexcept { return x > y; };
const auto binary_op_LE  = [](auto x, auto y) noexcept { return x <= y; };
const auto binary_op_GE  = [](auto x, auto y) noexcept { return x >= y; };
const auto binary_op_NE  = [](auto x, auto y) noexcept { return is_not_equal(x, y); };
const auto binary_op_EQ  = [](auto x, auto y) noexcept { return is_equal(x, y); };
const auto binary_op_LEG = [](auto x, auto y) noexcept { return (x < y) ? -1.0 : (x > y); };
const auto binary_op_AND = [](auto x, auto y) noexcept { return is_not_equal(x, 0) && is_not_equal(y, 0); };
const auto binary_op_OR  = [](auto x, auto y) noexcept { return is_not_equal(x, 0) || is_not_equal(y, 0); };
const auto binary_op_POW = [](double x, double y) noexcept { return std::pow(x, y); };
const auto binary_op_ADD = [](double x, double y) noexcept { return x + y; };
const auto binary_op_SUB = [](double x, double y) noexcept { return x - y; };
const auto binary_op_MUL = [](double x, double y) noexcept { return x * y; };
const auto binary_op_DIV = [](double x, double y) noexcept { return x / y; };
// clang-format on

//#define CHECK_UNUSED_VECTOR 1

#ifdef CHECK_UNUSED_VECTOR
// clang-format off
template <typename T>
class
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
CheckVector
{
 public:
  T *ptr;
  size_t m_count;

  CheckVector() : ptr(dummy), m_count(0) { }
  explicit CheckVector(size_t count) : ptr(dummy), m_count(count) { }
  CheckVector(size_t count, const T& value) : ptr(dummy), m_count(count) { ptr[0] = value; }

  T * begin() noexcept { return &ptr[0]; }
  T * end() noexcept { return &ptr[0] + 1; }
  const T * begin() const noexcept { return &ptr[0]; }
  const T * end() const noexcept { return &ptr[0] + 1; }

  bool empty() const { return true; }
  size_t size() const { return m_count; }
  void clear() { }
  void shrink_to_fit() { }

  void resize(size_t count) { m_count = count; }
  void resize(size_t count, const T& value) { m_count = count; ptr[0] = value; }

  void reserve(size_t new_cap) { (void)new_cap; }
  void push_back(const T& value) { ptr[0] = value; }
  void assign(size_t count, const T& value) { (void)count; ptr[0] = value; }

  T * data() noexcept { return ptr; }
  const T * data() const noexcept { return ptr; }

  T& operator[](size_t pos) { (void)pos; return ptr[0]; }
  const T & operator[](size_t pos) const { (void)pos; return ptr[0]; }
  CheckVector& operator=(const CheckVector& other) { (void)other; ptr = dummy; m_count = 0; return *this; }
  CheckVector& operator=(CheckVector&& other) { *this = other; ptr = dummy; m_count = 0; return *this; }
  CheckVector& operator=(const std::vector<T>& other) { (void)other; ptr = dummy; m_count = 0; return *this; }
  CheckVector& operator=(std::vector<T>&& other) { (void)other; ptr = dummy; m_count = 0; return *this; }
  // Copy constructor
  CheckVector(const CheckVector &v2) { m_count = 0; ptr = dummy; ptr[0] = v2.ptr[0]; }
  explicit CheckVector(std::vector<T> const &v2) { (void)v2; m_count = 0; ptr = dummy; }

  bool operator==(const CheckVector& other) const { (void)other; return true; }

 private:
  T dummy[1] {};
};
// clang-format on

template <typename T>
using Varray = CheckVector<T>;

#else

template <typename T>
using Varray = std::vector<T>;

#endif

template <typename T>
using Varray2D = Varray<Varray<T>>;

template <typename T>
using Varray3D = Varray2D<Varray<T>>;

template <typename T>
using Varray4D = Varray3D<Varray<T>>;

using Vmask = Varray<int8_t>;

// clang-format off
struct MinMax
{
  double min {std::numeric_limits<double>::max()};
  double max {-std::numeric_limits<double>::max()};
  size_t n {0};
  MinMax() {};
  MinMax(double rmin, double rmax, size_t rn) : min(rmin), max(rmax), n(rn) {};
  MinMax(double rmin, double rmax) : min(rmin), max(rmax), n(0) {};
};

struct MinMaxSum : MinMax
{
  double sum {0.0};
  MinMaxSum() {};
  MinMaxSum(double rmin, double rmax, double rsum, size_t rn) : sum(rsum) { min = rmin; max = rmax; n = rn; };
  MinMaxSum(double rmin, double rmax, double rsum) : sum(rsum) { min = rmin; max = rmax; n = 0; };
};

struct MinMaxMean : MinMax
{
  double mean {0.0};
  MinMaxMean() {};
  MinMaxMean(double rmin, double rmax, double rmean, size_t rn) : mean(rmean) { min = rmin; max = rmax; n = rn; };
  MinMaxMean(double rmin, double rmax, double rmean) : mean(rmean) { min = rmin; max = rmax; n = 0; };
};
// clang-format on

template <typename T>
inline void
varray_free(Varray<T> &v)
{
  v.clear();
  v.shrink_to_fit();
}

#define varrayResize(p, s) varray_resize((p), (s), __FILE__, __LINE__)
#define varrayResizeInit(p, s, v) varray_resize((p), (s), (v), __FILE__, __LINE__)

template <typename T>
inline void
varray_resize(Varray<T> &v, size_t count, const char *file, int line)
{
  if (v.size() != count)
  {
    try
    {
      v.resize(count);
    }
    catch (const std::exception &e)
    {
      std::cerr << "Exception caught when trying to allocate " << count << " vector elements: " << e.what() << " in " << file << ":"
                << line << '\n';
      throw;
    }
  }
}

template <typename T>
inline void
varray_resize(Varray<T> &v, size_t count, T value, const char *file, int line)
{
  try
  {
    v.resize(count, value);
  }
  catch (const std::exception &e)
  {
    std::cerr << "Exception caught when trying to allocate " << count << " vector elements: " << e.what() << " in " << file << ":"
              << line << '\n';
    throw;
  }
}

/*
template <class T, size_t ROW, size_t COL>
using Matrix = std::array<std::array<T, COL>, ROW>;
*/

template <typename T>
MinMaxSum varray_min_max_sum(Varray<T> const &v, size_t n, const MinMaxSum &mms);

template <typename T>
MinMaxSum varray_min_max_sum_mv(Varray<T> const &v, size_t n, const MinMaxSum &mms, double missval);

template <typename T>
MinMaxMean varray_min_max_mean(Varray<T> const &v, size_t n);

template <typename T>
MinMaxMean varray_min_max_mean_mv(Varray<T> const &v, size_t n, T missval);

template <typename T>
MinMax array_min_max_mask(const T *array, size_t n, Vmask const &mask);

void array_add_array(size_t n, double *array1, const double *array2);
void array_add_array_mv(size_t n, double *array1, const double *array2, double missval);

template <typename T1, typename T2>
void
array_copy(size_t n, const T1 *array1, T2 *array2)
{
  for (size_t i = 0; i < n; ++i) { array2[i] = array1[i]; }
}

template <typename T1, typename T2>
void
varray_copy(size_t n, Varray<T1> const &v1, Varray<T2> &v2)
{
  for (size_t i = 0; i < n; ++i) { v2[i] = v1[i]; }
}

template <typename T1, typename T2>
void
varray_divc(size_t n, Varray<T1> &v, T2 value)
{
  T1 c = value;
  for (size_t i = 0; i < n; ++i) { v[i] /= c; }
}

template <typename T, class UnaryOperation>
void
varray_transform(Varray<T> const &vIn, Varray<T> &vOut, UnaryOperation unary_op)
{
  assert(vIn.size() > 0);
  assert(vOut.size() > 0);
  assert(vOut.size() <= vIn.size());

  auto n = vIn.size();
  for (size_t i = 0; i < n; ++i) { vOut[i] = unary_op(vIn[i]); }
}

template <typename T>
size_t array_num_mv(size_t n, const T *array, double missval);

template <typename T>
size_t varray_num_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
MinMax varray_min_max(size_t n, Varray<T> const &v);

template <typename T>
MinMax varray_min_max(size_t n, const T *array);

template <typename T>
MinMax varray_min_max(Varray<T> const &v);

template <typename T>
MinMax varray_min_max_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
MinMax varray_min_max_mv(size_t n, const T *array, double missval);

template <typename T>
T varray_min(size_t n, Varray<T> const &v);

template <typename T>
T varray_max(size_t n, Varray<T> const &v);

template <typename T>
T varray_range(size_t n, Varray<T> const &v);

template <typename T>
double varray_min_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
double varray_max_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
double varray_range_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
double varray_sum(size_t n, Varray<T> const &v);

template <typename T>
double varray_sum_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
double varray_mean(size_t n, Varray<T> const &v);

template <typename T>
double varray_mean_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
double varray_weighted_mean(size_t n, Varray<T> const &v, Varray<double> const &w, double missval);

template <typename T>
double varray_weighted_mean_mv(size_t n, Varray<T> const &v, Varray<double> const &w, double missval);

template <typename T>
double varray_avg_mv(size_t n, Varray<T> const &v, double missval);

template <typename T>
double varray_weighted_avg_mv(size_t n, Varray<T> const &v, Varray<double> const &w, double missval);

template <typename T>
double varray_var(size_t n, Varray<T> const &v, size_t numMissVals, double missval);

template <typename T>
double varray_var_1(size_t n, Varray<T> const &v, size_t numMissVals, double missval);

template <typename T>
double varray_weighted_var(size_t n, Varray<T> const &v, Varray<double> const &w, size_t numMissVals, double missval);

template <typename T>
double varray_weighted_var_1(size_t n, Varray<T> const &v, Varray<double> const &w, size_t numMissVals, double missval);

template <typename T>
double varray_skew(size_t n, Varray<T> const &v, size_t numMissVals, double missval);

template <typename T>
double varray_kurt(size_t n, Varray<T> const &v, size_t numMissVals, double missval);

template <typename T>
double varray_median(size_t n, Varray<T> const &v, size_t numMissVals, double missval);

template <typename T>
double varray_count(size_t n, Varray<T> const &v, size_t numMissVals, double missval);

#endif  // VARRAY_H
