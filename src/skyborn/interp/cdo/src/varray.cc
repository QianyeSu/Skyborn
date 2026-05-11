/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cassert>

#include "compare.h"
#include "varray.h"
#include "cdo_omp.h"
#include "arithmetic.h"

template <typename T>
inline T
min_value(T v1, T v2)
{
  return (v1 < v2) ? v1 : v2;
}

template <typename T>
inline T
max_value(T v1, T v2)
{
  return (v1 > v2) ? v1 : v2;
}

template <typename T>
MinMax
varray_min_max_mv(size_t n, const T *array, double mv)
{
  T missval = mv;

  auto f_minmax_mv = [](auto a, auto mv_a, auto &vmin, auto &vmax, auto &nvals, auto is_NE)
  {
    if (is_NE(a, mv_a))
    {
      vmin = min_value(vmin, a);
      vmax = max_value(vmax, a);
      nvals++;
    }
  };

  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

  size_t nvals = 0;
  if (std::isnan(missval))
    for (size_t i = 0; i < n; ++i) { f_minmax_mv(array[i], missval, vmin, vmax, nvals, fp_is_not_equal); }
  else
    for (size_t i = 0; i < n; ++i) { f_minmax_mv(array[i], missval, vmin, vmax, nvals, is_not_equal); }

  return MinMax(vmin, vmax, nvals);
}

// Explicit instantiation
template MinMax varray_min_max_mv(size_t n, const float *array, double missval);
template MinMax varray_min_max_mv(size_t n, const double *array, double missval);

template <typename T>
MinMax
varray_min_max_mv(size_t n, Varray<T> const &v, double missval)
{
  return varray_min_max_mv(n, v.data(), missval);
}

// Explicit instantiation
template MinMax varray_min_max_mv(size_t n, Varray<float> const &v, double missval);
template MinMax varray_min_max_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
MinMaxSum
varray_min_max_sum(Varray<T> const &v, size_t n, const MinMaxSum &mms)
{
  auto f_minmaxsum = [](auto val, auto &vmin, auto &vmax, auto &vsum)
  {
    vmin = min_value(vmin, val);
    vmax = max_value(vmax, val);
    vsum += val;
  };

  auto vmin = mms.min;
  auto vmax = mms.max;
  auto vsum = mms.sum;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax) reduction(+ : vsum)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_minmaxsum((double) v[i], vmin, vmax, vsum); }
  }
  else
  {
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(min : vmin) reduction(max : vmax) reduction(+ : vsum)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_minmaxsum((double) v[i], vmin, vmax, vsum); }
  }

  return MinMaxSum(vmin, vmax, vsum, n);
}

// Explicit instantiation
template MinMaxSum varray_min_max_sum(Varray<float> const &v, size_t n, const MinMaxSum &mms);
template MinMaxSum varray_min_max_sum(Varray<double> const &v, size_t n, const MinMaxSum &mms);

template <typename T>
MinMaxSum
varray_min_max_sum_mv(Varray<T> const &v, size_t n, const MinMaxSum &mms, double missval)
{
  auto f_minmaxsum_mv = [](auto val, auto mv, auto &vmin, auto &vmax, auto &vsum, auto &nvals, auto is_NE)
  {
    if (is_NE(val, mv))
    {
      vmin = min_value(vmin, val);
      vmax = max_value(vmax, val);
      vsum += val;
      nvals++;
    }
  };

  T mv = static_cast<T>(missval);
  auto vmin = mms.min;
  auto vmax = mms.max;
  auto vsum = mms.sum;

  size_t nvals = 0;
  if (std::isnan(mv))
  {
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax) \
    reduction(+ : vsum, nvals)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_minmaxsum_mv((double) v[i], mv, vmin, vmax, vsum, nvals, fp_is_not_equal); }
  }
  else
  {
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax) \
    reduction(+ : vsum, nvals)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_minmaxsum_mv((double) v[i], mv, vmin, vmax, vsum, nvals, is_not_equal); }
  }

  if (nvals == 0 && is_equal(vmin, std::numeric_limits<double>::max())) vmin = mv;
  if (nvals == 0 && is_equal(vmax, -std::numeric_limits<double>::max())) vmax = mv;

  return MinMaxSum(vmin, vmax, vsum, nvals);
}

// Explicit instantiation
template MinMaxSum varray_min_max_sum_mv(Varray<float> const &v, size_t n, const MinMaxSum &mms, double missval);
template MinMaxSum varray_min_max_sum_mv(Varray<double> const &v, size_t n, const MinMaxSum &mms, double missval);

template <typename T>
MinMaxMean
varray_min_max_mean(Varray<T> const &v, size_t n)
{
  auto mms = varray_min_max_sum(v, n, MinMaxSum());
  auto rmean = (n != 0) ? mms.sum / static_cast<double>(n) : 0.0;
  return MinMaxMean(mms.min, mms.max, rmean, n);
}

// Explicit instantiation
template MinMaxMean varray_min_max_mean(Varray<float> const &v, size_t n);
template MinMaxMean varray_min_max_mean(Varray<double> const &v, size_t n);

template <typename T>
MinMaxMean
varray_min_max_mean_mv(Varray<T> const &v, size_t n, T missval)
{
  auto mms = varray_min_max_sum_mv(v, n, MinMaxSum(), missval);
  auto rmean = (mms.n != 0) ? mms.sum / static_cast<double>(mms.n) : missval;
  return MinMaxMean(mms.min, mms.max, rmean, mms.n);
}

// Explicit instantiation
template MinMaxMean varray_min_max_mean_mv(Varray<float> const &v, size_t n, float missval);
template MinMaxMean varray_min_max_mean_mv(Varray<double> const &v, size_t n, double missval);

template <typename T>
MinMax
array_min_max_mask(const T *const array, size_t n, Vmask const &mask)
{
  T rmin = std::numeric_limits<T>::max();
  T rmax = -std::numeric_limits<T>::max();

  if (!mask.empty())
  {
    for (size_t i = 0; i < n; ++i)
    {
      if (mask[i] == 0)
      {
        rmin = min_value(rmin, array[i]);
        rmax = max_value(rmax, array[i]);
      }
    }
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      rmin = min_value(rmin, array[i]);
      rmax = max_value(rmax, array[i]);
    }
  }

  return MinMax(rmin, rmax);
}

// Explicit instantiation
template MinMax array_min_max_mask(const float *const array, size_t n, Vmask const &mask);
template MinMax array_min_max_mask(const double *const array, size_t n, Vmask const &mask);

void
array_add_array(size_t n, double *array1, const double *array2)
{
#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
  for (size_t i = 0; i < n; ++i) { array1[i] += array2[i]; }
}

void
array_add_array_mv(size_t n, double *array1, const double *array2, double missval)
{
  if (std::isnan(missval))
  {
    for (size_t i = 0; i < n; ++i)
    {
      if (fp_is_not_equal(array2[i], missval)) { array1[i] = fp_is_equal(array1[i], missval) ? array2[i] : array1[i] + array2[i]; }
    }
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      if (is_not_equal(array2[i], missval)) { array1[i] = is_equal(array1[i], missval) ? array2[i] : array1[i] + array2[i]; }
    }
  }
}

auto count_mv = [](auto val, auto mv, auto &num, auto is_EQ)
{
  if (is_EQ(val, mv)) { num++; }
};

template <typename T>
size_t
array_num_mv(size_t n, const T *array, double mv)
{
  T missval = mv;
  size_t numMissVals = 0;

  if (std::isnan(missval))
  {
    for (size_t i = 0; i < n; ++i) { count_mv(array[i], missval, numMissVals, fp_is_equal); }
  }
  else
  {
    for (size_t i = 0; i < n; ++i) { count_mv(array[i], missval, numMissVals, is_equal); }
  }

  return numMissVals;
}

// Explicit instantiation
template size_t array_num_mv(size_t n, const float *array, double missval);
template size_t array_num_mv(size_t n, const double *array, double missval);

template <typename T>
size_t
varray_num_mv(size_t n, Varray<T> const &v, double mv)
{
  T missval = mv;

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  size_t numMissVals = 0;

  if (std::isnan(missval))
  {
    for (size_t i = 0; i < n; ++i) { count_mv(v[i], missval, numMissVals, fp_is_equal); }
  }
  else
  {
    for (size_t i = 0; i < n; ++i) { count_mv(v[i], missval, numMissVals, is_equal); }
  }

  return numMissVals;
}

// Explicit instantiation
template size_t varray_num_mv(size_t n, Varray<float> const &v, double missval);
template size_t varray_num_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
MinMax
varray_min_max(size_t n, const T *array)
{
  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
  for (size_t i = 0; i < n; ++i)
  {
    vmin = min_value(vmin, array[i]);
    vmax = max_value(vmax, array[i]);
  }

  return MinMax(vmin, vmax);
}

// Explicit instantiation
template MinMax varray_min_max(size_t n, const float *array);
template MinMax varray_min_max(size_t n, const double *array);

template <typename T>
MinMax
varray_min_max(size_t n, Varray<T> const &v)
{
  return varray_min_max(n, v.data());
}

// Explicit instantiation
template MinMax varray_min_max(size_t n, Varray<float> const &v);
template MinMax varray_min_max(size_t n, Varray<double> const &v);

template <typename T>
MinMax
varray_min_max(Varray<T> const &v)
{
  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

  auto n = v.size();
#ifndef __ICC  // wrong result with icc19
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (n > cdoMinLoopSize) default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
  for (size_t i = 0; i < n; ++i)
  {
    vmin = min_value(vmin, v[i]);
    vmax = max_value(vmax, v[i]);
  }

  return MinMax(vmin, vmax);
}

// Explicit instantiation
template MinMax varray_min_max(Varray<float> const &v);
template MinMax varray_min_max(Varray<double> const &v);

template <typename T>
T
varray_min(size_t n, Varray<T> const &v)
{
  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  auto vmin = v[0];
  if (n > cdoMinLoopSize)
  {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(min : vmin)
#endif
    for (size_t i = 0; i < n; ++i) { vmin = min_value(vmin, v[i]); }
  }
  else
  {
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(min : vmin)
#endif
    for (size_t i = 0; i < n; ++i) { vmin = min_value(vmin, v[i]); }
  }

  return vmin;
}

// Explicit instantiation
template float varray_min(size_t n, Varray<float> const &v);
template double varray_min(size_t n, Varray<double> const &v);

template <typename T>
T
varray_max(size_t n, Varray<T> const &v)
{
  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  auto vmax = v[0];
  if (n > cdoMinLoopSize)
  {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(max : vmax)
#endif
    for (size_t i = 0; i < n; ++i) { vmax = max_value(vmax, v[i]); }
  }
  else
  {
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(max : vmax)
#endif
    for (size_t i = 0; i < n; ++i) { vmax = max_value(vmax, v[i]); }
  }

  return vmax;
}

// Explicit instantiation
template float varray_max(size_t n, Varray<float> const &v);
template double varray_max(size_t n, Varray<double> const &v);

template <typename T>
T
varray_range(size_t n, Varray<T> const &v)
{
  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  auto vmin = v[0];
  auto vmax = v[0];
  if (n > cdoMinLoopSize)
  {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
    for (size_t i = 0; i < n; ++i)
    {
      vmin = min_value(vmin, v[i]);
      vmax = max_value(vmax, v[i]);
    }
  }
  else
  {
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(min : vmin) reduction(max : vmax)
#endif
    for (size_t i = 0; i < n; ++i)
    {
      vmin = min_value(vmin, v[i]);
      vmax = max_value(vmax, v[i]);
    }
  }

  return (vmax - vmin);
}

// Explicit instantiation
template float varray_range(size_t n, Varray<float> const &v);
template double varray_range(size_t n, Varray<double> const &v);

template <typename T>
double
varray_min_mv(size_t n, Varray<T> const &v, double mv)
{
  T missval = mv;

  auto f_min_mv = [](auto a, auto mv_a, auto &vmin, auto is_NE)
  {
    if (is_NE(a, mv_a)) { vmin = min_value(vmin, a); }
  };

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  T vmin = std::numeric_limits<T>::max();

  if (std::isnan(missval))
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(min : vmin)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_min_mv(v[i], missval, vmin, fp_is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(min : vmin)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_min_mv(v[i], missval, vmin, fp_is_not_equal); }
    }
  }
  else
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(min : vmin)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_min_mv(v[i], missval, vmin, is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(min : vmin)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_min_mv(v[i], missval, vmin, is_not_equal); }
    }
  }

  if (is_equal(vmin, std::numeric_limits<T>::max())) vmin = missval;

  return vmin;
}

// Explicit instantiation
template double varray_min_mv(size_t n, Varray<float> const &v, double missval);
template double varray_min_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
double
varray_max_mv(size_t n, Varray<T> const &v, double mv)
{
  T missval = mv;

  auto f_max_mv = [](auto a, auto mv_a, auto &vmax, auto is_NE)
  {
    if (is_NE(a, mv_a)) { vmax = max_value(vmax, a); }
  };

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  T vmax = -std::numeric_limits<T>::max();

  if (std::isnan(missval))
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_max_mv(v[i], missval, vmax, fp_is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_max_mv(v[i], missval, vmax, fp_is_not_equal); }
    }
  }
  else
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_max_mv(v[i], missval, vmax, is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_max_mv(v[i], missval, vmax, is_not_equal); }
    }
  }

  if (is_equal(vmax, -std::numeric_limits<T>::max())) vmax = missval;

  return vmax;
}

// Explicit instantiation
template double varray_max_mv(size_t n, Varray<float> const &v, double missval);
template double varray_max_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
double
varray_range_mv(size_t n, Varray<T> const &v, double mv)
{
  T missval = mv;

  auto f_minmax_mv = [](auto a, auto mv_a, auto &vmin, auto &vmax, auto is_NE)
  {
    if (is_NE(a, mv_a))
    {
      vmin = min_value(vmin, a);
      vmax = max_value(vmax, a);
    }
  };

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  T vmin = std::numeric_limits<T>::max();
  T vmax = -std::numeric_limits<T>::max();

  if (std::isnan(missval))
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_minmax_mv(v[i], missval, vmin, vmax, fp_is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(min : vmin) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_minmax_mv(v[i], missval, vmin, vmax, fp_is_not_equal); }
    }
  }
  else
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(min : vmin) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_minmax_mv(v[i], missval, vmin, vmax, is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(min : vmin) reduction(max : vmax)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_minmax_mv(v[i], missval, vmin, vmax, is_not_equal); }
    }
  }

  return (is_equal(vmin, std::numeric_limits<T>::max()) && is_equal(vmax, -std::numeric_limits<T>::max())) ? missval : vmax - vmin;
}

// Explicit instantiation
template double varray_range_mv(size_t n, Varray<float> const &v, double missval);
template double varray_range_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
double
varray_sum(size_t n, Varray<T> const &v)
{
  // assert(n > 0); // failed in remapcon
  assert(v.size() > 0);
  assert(n <= v.size());

  double sum = 0.0;
  if (n > cdoMinLoopSize)
  {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : sum)
#endif
    for (size_t i = 0; i < n; ++i) { sum += v[i]; }
  }
  else
  {
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sum)
#endif
    for (size_t i = 0; i < n; ++i) { sum += v[i]; }
  }

  return sum;
}

// Explicit instantiation
template double varray_sum(size_t n, Varray<float> const &v);
template double varray_sum(size_t n, Varray<double> const &v);

template <typename T>
double
varray_sum_mv(size_t n, Varray<T> const &v, double mv)
{
  T missval = mv;

  auto f_sum_mv = [](auto a, auto mv_a, auto &sum, auto &nvals, auto is_NE)
  {
    if (is_NE(a, mv_a))
    {
      sum += a;
      nvals++;
    }
  };

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  double sum = 0.0;
  size_t nvals = 0;

  if (std::isnan(missval))
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : sum, nvals)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_sum_mv(v[i], missval, sum, nvals, fp_is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sum, nvals)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_sum_mv(v[i], missval, sum, nvals, fp_is_not_equal); }
    }
  }
  else
  {
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : sum, nvals)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_sum_mv(v[i], missval, sum, nvals, is_not_equal); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sum, nvals)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_sum_mv(v[i], missval, sum, nvals, is_not_equal); }
    }
  }

  if (!nvals) { sum = missval; }

  return sum;
}

// Explicit instantiation
template double varray_sum_mv(size_t n, Varray<float> const &v, double missval);
template double varray_sum_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
double
varray_mean(size_t n, Varray<T> const &v)
{
  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  auto sum = varray_sum(n, v);

  return sum / n;
}

// Explicit instantiation
template double varray_mean(size_t n, Varray<float> const &v);
template double varray_mean(size_t n, Varray<double> const &v);

template <typename T>
double
varray_mean_mv(size_t n, Varray<T> const &v, double mv)
{
  T missval = mv;

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  auto is_NE = fp_is_not_equal;
  auto is_EQ = fp_is_equal;
  double sum = 0.0, sumw = 0.0;

  for (size_t i = 0; i < n; ++i)
    if (is_NE(v[i], missval))
    {
      sum += v[i];
      sumw += 1;
    }

  double missval1 = missval, missval2 = missval;
  return DIVM(sum, sumw);
}

// Explicit instantiation
template double varray_mean_mv(size_t n, Varray<float> const &v, double missval);
template double varray_mean_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
double
varray_weighted_mean(size_t n, Varray<T> const &v, Varray<double> const &w, double mv)
{
  T missval = mv;

  auto f_weighted_mean = [](auto aw, auto a, auto &sum, auto &sumw)
  {
    sum += aw * a;
    sumw += aw;
  };

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());
  assert(n <= w.size());

  double sum = 0.0, sumw = 0.0;
  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : sum, sumw)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_weighted_mean(w[i], v[i], sum, sumw); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sum, sumw)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_weighted_mean(w[i], v[i], sum, sumw); }
  }

  return is_equal(sumw, 0.0) ? missval : sum / sumw;
}

// Explicit instantiation
template double varray_weighted_mean(size_t n, Varray<float> const &v, Varray<double> const &w, double missval);
template double varray_weighted_mean(size_t n, Varray<double> const &v, Varray<double> const &w, double missval);

template <typename T>
double
varray_weighted_mean_mv(size_t n, Varray<T> const &v, Varray<double> const &w, double mv)
{
  T missval = mv;

  auto f_weighted_mean_mv = [](auto aw, auto a, auto mv_a, auto &sum, auto &sumw, auto is_NE)
  {
    if (is_NE(a, mv_a) && is_NE(aw, mv_a))
    {
      sum += aw * a;
      sumw += aw;
    }
  };

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());
  assert(n <= w.size());

  double missval1 = missval, missval2 = missval;
  double sum = 0.0, sumw = 0.0;
  double wmean = 0.0;

  if (std::isnan(missval))
  {
    auto is_NE = fp_is_not_equal;
    auto is_EQ = fp_is_equal;
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for default(shared) schedule(static) reduction(+ : sum, sumw)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_weighted_mean_mv(w[i], v[i], missval1, sum, sumw, is_NE); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sum, sumw)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_weighted_mean_mv(w[i], v[i], missval1, sum, sumw, is_NE); }
    }
    wmean = DIVM(sum, sumw);
  }
  else
  {
    auto is_NE = is_not_equal;
    auto is_EQ = is_equal;
    if (n > cdoMinLoopSize)
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : sum, sumw)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_weighted_mean_mv(w[i], v[i], missval1, sum, sumw, is_NE); }
    }
    else
    {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : sum, sumw)
#endif
#endif
      for (size_t i = 0; i < n; ++i) { f_weighted_mean_mv(w[i], v[i], missval1, sum, sumw, is_NE); }
    }
    wmean = DIVM(sum, sumw);
  }

  return wmean;
}

// Explicit instantiation
template double varray_weighted_mean_mv(size_t n, Varray<float> const &v, Varray<double> const &w, double missval);
template double varray_weighted_mean_mv(size_t n, Varray<double> const &v, Varray<double> const &w, double missval);

template <typename T>
double
varray_avg_mv(size_t n, Varray<T> const &v, double mv)
{
  T missval = mv;

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  auto is_EQ = fp_is_equal;
  double missval1 = missval, missval2 = missval;
  double sum = 0.0, sumw = 0.0;

  for (size_t i = 0; i < n; ++i)
  {
    sum = ADDM(sum, v[i]);
    sumw += 1;
  }

  return DIVM(sum, sumw);
}

// Explicit instantiation
template double varray_avg_mv(size_t n, Varray<float> const &v, double missval);
template double varray_avg_mv(size_t n, Varray<double> const &v, double missval);

template <typename T>
double
varray_weighted_avg_mv(size_t n, Varray<T> const &v, Varray<double> const &w, double mv)
{
  T missval = mv;

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());
  assert(n <= w.size());

  auto is_EQ = fp_is_equal;
  double missval1 = missval, missval2 = missval;
  double sum = 0.0, sumw = 0.0;

  for (size_t i = 0; i < n; ++i)
  {
    if (fp_is_not_equal(w[i], missval))
    {
      sum = ADDM(sum, MULM(w[i], v[i]));
      sumw = ADDM(sumw, w[i]);
    }
  }

  return DIVM(sum, sumw);
}

// Explicit instantiation
template double varray_weighted_avg_mv(size_t n, Varray<float> const &v, Varray<double> const &w, double missval);
template double varray_weighted_avg_mv(size_t n, Varray<double> const &v, Varray<double> const &w, double missval);

template <typename T>
static void
varray_prevarsum0(size_t n, Varray<T> const &v, double &rsum, double &rsumw)
{
  rsum = 0.0;
  if (n > cdoMinLoopSize)
  {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum)
#endif
    for (size_t i = 0; i < n; ++i) { rsum += v[i]; }
  }
  else
  {
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum)
#endif
    for (size_t i = 0; i < n; ++i) { rsum += v[i]; }
  }
  rsumw = static_cast<double>(n);
}

template <typename T>
static void
varray_prevarsum0_mv(size_t n, Varray<T> const &v, double missval, double &rsum, double &rsumw)
{
  auto f_prevarsum0_mv = [](auto a, auto mv_a, auto &sum, auto &sumw)
  {
    if (fp_is_not_equal(a, mv_a))
    {
      sum += a;
      sumw += 1.0;
    }
  };

  rsum = rsumw = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum, rsumw)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prevarsum0_mv(v[i], missval, rsum, rsumw); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum, rsumw)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prevarsum0_mv(v[i], missval, rsum, rsumw); }
  }
}

template <typename T>
static void
varray_prevarsum(size_t n, Varray<T> const &v, double &rsum, double &rsumw, double &rsumq, double &rsumwq)
{
  auto f_prevarsum = [](auto a, auto &sum, auto &sumq)
  {
    sum += a;
    sumq += a * a;
  };

  rsum = rsumq = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum, rsumq)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prevarsum((double) v[i], rsum, rsumq); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum, rsumq)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prevarsum((double) v[i], rsum, rsumq); }
  }

  rsumw = static_cast<double>(n);
  rsumwq = static_cast<double>(n);
}

template <typename T>
static void
varray_prevarsum_mv(size_t n, Varray<T> const &v, T missval, double &rsum, double &rsumw, double &rsumq, double &rsumwq)
{
  auto f_prevarsum = [](auto a, auto mv_a, auto &sum, auto &sumq, auto &sumw, auto &sumwq)
  {
    if (fp_is_not_equal(a, mv_a))
    {
      double ad = (double) a;
      sum += ad;
      sumq += ad * ad;
      sumw += 1.0;
      sumwq += 1.0;
    }
  };

  rsum = rsumq = rsumw = rsumwq = 0.0;

#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared) schedule(static) reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
  for (size_t i = 0; i < n; ++i) { f_prevarsum(v[i], missval, rsum, rsumq, rsumw, rsumwq); }
}

template <typename T>
double
varray_var(size_t n, Varray<T> const &v, size_t numMissVals, double mv)
{
  T missval = mv;

  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (numMissVals > 0) { varray_prevarsum_mv(n, v, missval, rsum, rsumw, rsumq, rsumwq); }
  else { varray_prevarsum(n, v, rsum, rsumw, rsumq, rsumwq); }

  auto rvar = is_not_equal(rsumw, 0.0) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) { rvar = 0.0; }

  return rvar;
}

// Explicit instantiation
template double varray_var(size_t n, Varray<float> const &v, size_t numMissVals, double missval);
template double varray_var(size_t n, Varray<double> const &v, size_t numMissVals, double missval);

template <typename T>
double
varray_var_1(size_t n, Varray<T> const &v, size_t numMissVals, double mv)
{
  T missval = mv;

  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (numMissVals > 0) { varray_prevarsum_mv(n, v, missval, rsum, rsumw, rsumq, rsumwq); }
  else { varray_prevarsum(n, v, rsum, rsumw, rsumq, rsumwq); }

  auto rvar = (rsumw * rsumw > rsumwq) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw - rsumwq) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) { rvar = 0.0; }

  return rvar;
}

// Explicit instantiation
template double varray_var_1(size_t n, Varray<float> const &v, size_t numMissVals, double missval);
template double varray_var_1(size_t n, Varray<double> const &v, size_t numMissVals, double missval);

template <typename T>
static void
varray_weighted_prevarsum(size_t n, Varray<T> const &v, Varray<double> const &w, double &rsum, double &rsumw, double &rsumq,
                          double &rsumwq)
{
  auto f_weighted_prevarsum = [](auto aw, auto a, auto &sum, auto &sumq, auto &sumw, auto &sumwq)
  {
    sum += aw * a;
    sumq += aw * a * a;
    sumw += aw;
    sumwq += aw * aw;
  };

  rsum = rsumq = rsumw = rsumwq = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_weighted_prevarsum(w[i], (double) v[i], rsum, rsumq, rsumw, rsumwq); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_weighted_prevarsum(w[i], (double) v[i], rsum, rsumq, rsumw, rsumwq); }
  }
}

template <typename T>
static void
varray_weighted_prevarsum_mv(size_t n, Varray<T> const &v, Varray<double> const &w, double missval, double &rsum, double &rsumw,
                             double &rsumq, double &rsumwq)
{
  auto f_weighted_prevarsum_mv = [](auto aw, auto a, auto mv_a, auto &sum, auto &sumq, auto &sumw, auto &sumwq)
  {
    if (fp_is_not_equal(a, mv_a) && fp_is_not_equal(aw, mv_a))
    {
      sum += aw * a;
      sumq += aw * a * a;
      sumw += aw;
      sumwq += aw * aw;
    }
  };

  rsum = rsumq = rsumw = rsumwq = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_weighted_prevarsum_mv(w[i], (double) v[i], missval, rsum, rsumq, rsumw, rsumwq); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum, rsumq, rsumw, rsumwq)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_weighted_prevarsum_mv(w[i], (double) v[i], missval, rsum, rsumq, rsumw, rsumwq); }
  }
}

template <typename T>
double
varray_weighted_var(size_t n, Varray<T> const &v, Varray<double> const &w, size_t numMissVals, double mv)
{
  T missval = mv;

  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (numMissVals > 0) { varray_weighted_prevarsum_mv(n, v, w, missval, rsum, rsumw, rsumq, rsumwq); }
  else { varray_weighted_prevarsum(n, v, w, rsum, rsumw, rsumq, rsumwq); }

  double rvar = is_not_equal(rsumw, 0) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) { rvar = 0.0; }

  return rvar;
}

// Explicit instantiation
template double varray_weighted_var(size_t n, Varray<float> const &v, Varray<double> const &w, size_t numMissVals, double missval);
template double varray_weighted_var(size_t n, Varray<double> const &v, Varray<double> const &w, size_t numMissVals, double missval);

template <typename T>
double
varray_weighted_var_1(size_t n, Varray<T> const &v, Varray<double> const &w, size_t numMissVals, double mv)
{
  T missval = mv;

  double rsum = 0.0, rsumw = 0.0, rsumq = 0.0, rsumwq = 0.0;
  if (numMissVals > 0) { varray_weighted_prevarsum_mv(n, v, w, missval, rsum, rsumw, rsumq, rsumwq); }
  else { varray_weighted_prevarsum(n, v, w, rsum, rsumw, rsumq, rsumwq); }

  double rvar = (rsumw * rsumw > rsumwq) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw - rsumwq) : missval;
  if (rvar < 0.0 && rvar > -1.e-5) { rvar = 0.0; }

  return rvar;
}

// Explicit instantiation
template double varray_weighted_var_1(size_t n, Varray<float> const &v, Varray<double> const &w, size_t numMissVals,
                                      double missval);
template double varray_weighted_var_1(size_t n, Varray<double> const &v, Varray<double> const &w, size_t numMissVals,
                                      double missval);

template <typename T>
static void
varray_prekurtsum(size_t n, Varray<T> const &v, double mean, double &rsum3w, double &rsum2diff, double &rsum4diff)
{
  auto f_prekurtsum = [](auto vdiff, auto &sum2diff, auto &sum4diff)
  {
    sum2diff += vdiff * vdiff;
    sum4diff += vdiff * vdiff * vdiff * vdiff;
  };

  rsum2diff = rsum4diff = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum2diff, rsum4diff)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prekurtsum(v[i] - mean, rsum2diff, rsum4diff); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum2diff, rsum4diff)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prekurtsum(v[i] - mean, rsum2diff, rsum4diff); }
  }

  rsum3w = static_cast<double>(n);
}

template <typename T>
static void
varray_prekurtsum_mv(size_t n, Varray<T> const &v, T missval, double mean, double &rsum3w, double &rsum2diff, double &rsum4diff)
{
  auto f_prekurtsum_mv = [](auto a, auto mv_a, auto meanval, auto &sum2diff, auto &sum4diff, auto &sum3w)
  {
    if (fp_is_not_equal(a, mv_a))
    {
      double vdiff = a - meanval;
      sum2diff += vdiff * vdiff;
      sum4diff += vdiff * vdiff * vdiff * vdiff;
      sum3w += 1;
    }
  };

  rsum3w = rsum2diff = rsum4diff = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum2diff, rsum4diff, rsum3w)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prekurtsum_mv(v[i], missval, mean, rsum2diff, rsum4diff, rsum3w); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum2diff, rsum4diff, rsum3w)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_prekurtsum_mv(v[i], missval, mean, rsum2diff, rsum4diff, rsum3w); }
  }
}

template <typename T>
double
varray_kurt(size_t n, Varray<T> const &v, size_t numMissVals, double mv)
{
  T missval = mv;

  double rsum3w{};  // 3rd moment variables
  double rsum2diff{}, rsum4diff{};
  double rsum{}, rsumw{};

  if (numMissVals > 0)
  {
    varray_prevarsum0_mv(n, v, missval, rsum, rsumw);
    varray_prekurtsum_mv(n, v, missval, (rsum / rsumw), rsum3w, rsum2diff, rsum4diff);
  }
  else
  {
    varray_prevarsum0(n, v, rsum, rsumw);
    varray_prekurtsum(n, v, (rsum / rsumw), rsum3w, rsum2diff, rsum4diff);
  }

  if (is_equal(rsum3w, 0.0) || is_equal(rsum2diff, 0.0)) { return missval; }

  auto rkurt = ((rsum4diff / rsum3w) / std::pow(rsum2diff / rsum3w, 2)) - 3.0;
  if (rkurt < 0.0 && rkurt > -1.e-5) { rkurt = 0.0; }

  return rkurt;
}

// Explicit instantiation
template double varray_kurt(size_t n, Varray<float> const &v, size_t numMissVals, double missval);
template double varray_kurt(size_t n, Varray<double> const &v, size_t numMissVals, double missval);

template <typename T>
static void
varray_preskewsum(size_t n, Varray<T> const &v, double mean, double &rsum3w, double &rsum3diff, double &rsum2diff)
{
  auto f_preskewsum = [](auto vdiff, auto &sum3diff, auto &sum2diff)
  {
    sum3diff += vdiff * vdiff * vdiff;
    sum2diff += vdiff * vdiff;
  };

  rsum2diff = 0.0;
  rsum3diff = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum2diff, rsum3diff)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_preskewsum(v[i] - mean, rsum3diff, rsum2diff); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum2diff, rsum3diff)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_preskewsum(v[i] - mean, rsum3diff, rsum2diff); }
  }

  rsum3w = static_cast<double>(n);
}

template <typename T>
static void
varray_preskewsum_mv(size_t n, Varray<T> const &v, T missval, double mean, double &rsum3w, double &rsum3diff, double &rsum2diff)
{
  auto f_preskewsum_mv = [](auto a, auto mv_a, auto meanval, auto &sum3diff, auto &sum2diff, auto &sum3w)
  {
    if (fp_is_not_equal(a, mv_a))
    {
      double vdiff = a - meanval;
      sum3diff += vdiff * vdiff * vdiff;
      sum2diff += vdiff * vdiff;
      sum3w += 1;
    }
  };

  rsum3w = rsum3diff = rsum2diff = 0.0;

  if (n > cdoMinLoopSize)
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(shared) schedule(static) reduction(+ : rsum2diff, rsum3diff, rsum3w)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_preskewsum_mv(v[i], missval, mean, rsum3diff, rsum2diff, rsum3w); }
  }
  else
  {
#ifndef __ICC  // internal error with icc22: lambda not supported
#ifdef HAVE_OPENMP4
#pragma omp simd reduction(+ : rsum2diff, rsum3diff, rsum3w)
#endif
#endif
    for (size_t i = 0; i < n; ++i) { f_preskewsum_mv(v[i], missval, mean, rsum3diff, rsum2diff, rsum3w); }
  }
}

template <typename T>
double
varray_skew(size_t n, Varray<T> const &v, size_t numMissVals, double mv)
{
  T missval = mv;

  double rsum3w{};  // 3rd moment variables
  double rsum3diff{}, rsum2diff{};
  double rsum{}, rsumw{};

  if (numMissVals > 0)
  {
    varray_prevarsum0_mv(n, v, missval, rsum, rsumw);
    varray_preskewsum_mv(n, v, missval, (rsum / rsumw), rsum3w, rsum3diff, rsum2diff);
  }
  else
  {
    varray_prevarsum0(n, v, rsum, rsumw);
    varray_preskewsum(n, v, (rsum / rsumw), rsum3w, rsum3diff, rsum2diff);
  }

  if (is_equal(rsum3w, 0.0) || is_equal(rsum3w, 1.0) || is_equal(rsum2diff, 0.0)) { return missval; }

  auto rskew = (rsum3diff / rsum3w) / std::pow((rsum2diff) / (rsum3w - 1.0), 1.5);
  if (rskew < 0.0 && rskew > -1.e-5) { rskew = 0.0; }

  return rskew;
}

// Explicit instantiation
template double varray_skew(size_t n, Varray<float> const &v, size_t numMissVals, double missval);
template double varray_skew(size_t n, Varray<double> const &v, size_t numMissVals, double missval);

#include <algorithm>

template <typename T>
static double
get_nth_element(T *array, size_t ngth, size_t n)
{
  std::nth_element(array, array + n, array + ngth);
  return array[n];
}

template <typename T>
static double
f_median(size_t n, Varray<T> &v)
{
  if (n % 2 == 0)
  {
    auto k = n / 2;
    auto vk1 = get_nth_element(v.data(), n, k - 1);
    auto vk2 = get_nth_element(v.data(), n, k);
    return (vk1 + vk2) * 0.5;
  }
  else
  {
    auto k = (n + 1) / 2;
    return get_nth_element(v.data(), n, k - 1);
  }
}

template <typename T>
double
varray_median(size_t n, Varray<T> const &v, size_t numMissVals, double mv)
{
  T missval = mv;

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  double median = missval;

  if (numMissVals == 0)
  {
    Varray<T> v2 = v;
    median = f_median(n, v2);
  }
  else
  {
    Varray<T> v2(n);
    size_t count = 0;
    for (size_t i = 0; i < n; ++i)
    {
      if (fp_is_not_equal(v[i], missval)) { v2[count++] = v[i]; }
    }
    if (count > 0 && count < n) { median = f_median(count, v2); }
  }

  return median;
}

// Explicit instantiation
template double varray_median(size_t n, Varray<float> const &v, size_t numMissVals, double missval);
template double varray_median(size_t n, Varray<double> const &v, size_t numMissVals, double missval);

template <typename T>
double
varray_count(size_t n, Varray<T> const &v, size_t numMissVals, double mv)
{
  T missval = mv;

  assert(n > 0);
  assert(v.size() > 0);
  assert(n <= v.size());

  size_t count = n;

  if (numMissVals > 0)
  {
    count = 0;
    for (size_t i = 0; i < n; ++i)
    {
      if (fp_is_not_equal(v[i], missval)) { count++; }
    }
  }

  return static_cast<double>(count);
}

// Explicit instantiation
template double varray_count(size_t n, Varray<float> const &v, size_t numMissVals, double missval);
template double varray_count(size_t n, Varray<double> const &v, size_t numMissVals, double missval);
