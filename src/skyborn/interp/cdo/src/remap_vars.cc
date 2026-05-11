/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <utility>

#include "field.h"
#include "cdo_timer.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_omp.h"
#include "remap_vars.h"

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere

  -----------------------------------------------------------------------
*/

template <typename T1, typename T2>
static void
remap_first_order(Varray<T2> &tgtArray, RemapVars const &rv, Varray<T1> const &srcArray)
{
  auto numLinks = rv.numLinks;
  auto numWeights = rv.numWeights;
  auto const &weights = rv.weights;
  auto const &tgtIndices = rv.tgtCellIndices;
  auto const &srcIndices = rv.srcCellIndices;
  auto numLinksPerValue = rv.numLinksPerValue;

  if (numLinksPerValue > 0 && numWeights == 1)
  {
    size_t nlinks = numLinks / numLinksPerValue;

    if (numLinksPerValue == 1)
    {
#ifdef _OPENMP
#pragma omp parallel for if (nlinks > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < nlinks; ++i) { tgtArray[tgtIndices[i]] = srcArray[srcIndices[i]] * weights[i]; }
    }
    else if (numLinksPerValue == 2)
    {
#ifdef _OPENMP
#pragma omp parallel for if (nlinks > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < nlinks; ++i)
      {
        auto noff = i * 2;
        const auto *const s = &srcIndices[noff];
        const auto *const w = &weights[noff];
        auto tgtPoint = srcArray[s[0]] * w[0] + srcArray[s[1]] * w[1];
        tgtArray[tgtIndices[noff]] = tgtPoint;
      }
    }
    else if (numLinksPerValue == 3)
    {
#ifdef _OPENMP
#pragma omp parallel for if (nlinks > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < nlinks; ++i)
      {
        auto noff = i * 3;
        const auto *const s = &srcIndices[noff];
        const auto *const w = &weights[noff];
        auto tgtPoint = srcArray[s[0]] * w[0] + srcArray[s[1]] * w[1] + srcArray[s[2]] * w[2];
        tgtArray[tgtIndices[noff]] = tgtPoint;
      }
    }
    else if (numLinksPerValue == 4)
    {
#ifdef _OPENMP
#pragma omp parallel for if (nlinks > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < nlinks; ++i)
      {
        auto noff = i * 4;
        const auto *const s = &srcIndices[noff];
        const auto *const w = &weights[noff];
        auto tgtPoint = srcArray[s[0]] * w[0] + srcArray[s[1]] * w[1] + srcArray[s[2]] * w[2] + srcArray[s[3]] * w[3];
        tgtArray[tgtIndices[noff]] = tgtPoint;
      }
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for if (nlinks > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < nlinks; ++i)
      {
        auto noff = i * numLinksPerValue;
        const auto *const s = &srcIndices[noff];
        const auto *const w = &weights[noff];
        auto tgtPoint = srcArray[s[0]] * w[0];
        for (size_t k = 1; k < (size_t) numLinksPerValue; ++k) { tgtPoint += srcArray[s[k]] * w[k]; }
        tgtArray[tgtIndices[noff]] = tgtPoint;
      }
    }
  }
  else
  {
    if (numWeights == 1 && rv.linksOffset.size() > 0 && rv.linksPerValue.size() > 0)
    {
      auto const &linksOffset = rv.linksOffset;
      auto const &linksPerValue = rv.linksPerValue;
      auto tgtGridSize = tgtArray.size();
#ifdef _OPENMP
#pragma omp parallel for if (tgtGridSize > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < tgtGridSize; ++i)
      {
        if (linksPerValue[i] > 0)
        {
          auto offset = linksOffset[i];
          auto nlinks = linksPerValue[i];
          tgtArray[i] = srcArray[srcIndices[offset]] * weights[offset];
          for (size_t k = 1; k < nlinks; ++k) { tgtArray[i] += srcArray[srcIndices[offset + k]] * weights[offset + k]; }
        }
      }
    }
    else
    {
      for (size_t i = 0; i < numLinks; ++i) { tgtArray[tgtIndices[i]] = static_cast<T2>(0.0); }
      for (size_t i = 0; i < numLinks; ++i) { tgtArray[tgtIndices[i]] += srcArray[srcIndices[i]] * weights[numWeights * i]; }
    }
  }
}

template <typename T1, typename T2>
static void
remap_second_order(Varray<T2> &tgtArray, RemapVars const &rv, Varray<T1> const &srcArray, RemapGradients &gradients)
{
  auto const &grad1 = gradients.lat;
  auto const &grad2 = gradients.lon;
  auto const &grad3 = gradients.latLon;

  auto numLinks = rv.numLinks;
  auto numWeights = rv.numWeights;
  auto const &weights = rv.weights;
  auto const &tgtIndices = rv.tgtCellIndices;
  auto const &srcIndices = rv.srcCellIndices;
  auto numLinksPerValue = rv.numLinksPerValue;

  if (numWeights == 3)
  {
    for (size_t i = 0; i < numLinks; ++i) { tgtArray[tgtIndices[i]] = static_cast<T2>(0.0); }
    for (size_t i = 0; i < numLinks; ++i)
    {
      auto k = srcIndices[i];
      const auto *const w = &weights[3 * i];
      tgtArray[tgtIndices[i]] += srcArray[k] * w[0] + grad1[k] * w[1] + grad2[k] * w[2];
      // printf("%zu %zu %.5f %.5f %.5f %.5f %.5f\n", i, k, grad1[k], grad2[k], w[0], w[1], w[2]);
    }
  }
  else if (numWeights == 4)
  {
    if (numLinksPerValue == 4)
    {
      size_t nlinks = numLinks / numLinksPerValue;
#ifdef _OPENMP
#pragma omp parallel for if (nlinks > cdoMinLoopSize) default(shared) schedule(static)
#endif
      for (size_t i = 0; i < nlinks; ++i)
      {
        double tgtPoint = 0.0;
        for (int k = 0; k < 4; ++k)
        {
          auto noff = i * 4 + k;
          auto ik = srcIndices[noff];
          const auto *const w = &weights[numWeights * noff];
          tgtPoint += srcArray[ik] * w[0] + grad1[ik] * w[1] + grad2[ik] * w[2] + grad3[ik] * w[3];
        }
        tgtArray[tgtIndices[i * 4]] = tgtPoint;
      }
    }
    else
    {
      for (size_t i = 0; i < numLinks; ++i) { tgtArray[tgtIndices[i]] = static_cast<T2>(0.0); }
      for (size_t i = 0; i < numLinks; ++i)
      {
        auto k = srcIndices[i];
        const auto *const w = &weights[4 * i];
        tgtArray[tgtIndices[i]] += srcArray[k] * w[0] + grad1[k] * w[1] + grad2[k] * w[2] + grad3[k] * w[3];
      }
    }
  }
}

template <typename T1, typename T2>
static void
remap(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double tgtMissval, size_t tgtSize, RemapVars const &rv,
      RemapGradients &gradients)
{
  T2 missval = tgtMissval;
  /*
    Input arrays:

      tgtIndices    target indices for each link
      srcIndices    source indices for each link
      numWeights    num of weights used in remapping
      weights       remapping weights for each link
      srcArray      array with source field to be remapped

    Optional:

      gradients  gradient arrays on source grid necessary for higher-order remappings

    Output variables:

      tgtArray  array for remapped field on target grid
  */
  if (Options::cdoVerbose) cdo_print("Remap links per value: %ld", rv.numLinksPerValue);

  cdo::timer timer;

  // Check the order of the interpolation

  auto firstOrder = (gradients.lat.size() == 0);

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd if (tgtSize > cdoMinLoopSize) default(shared) schedule(static)
#endif
  for (size_t n = 0; n < tgtSize; ++n) { tgtArray[n] = missval; }

  if (firstOrder)  // First order remapping
  {
    remap_first_order(tgtArray, rv, srcArray);
  }
  else  // Second order remapping
  {
    remap_second_order(tgtArray, rv, srcArray, gradients);
  }

  if (Options::cdoVerbose) cdo_print("Remap: %.2f seconds", timer.elapsed());
}

void
remap_field(Field &field2, double missval, size_t gridsize2, RemapVars const &rv, Field const &field1, RemapGradients &gradients)
{
  auto func = [&](auto const &v1, auto &v2) { remap(v1, v2, missval, gridsize2, rv, gradients); };
  field_operation2(func, field1, field2);
}

static size_t
get_max_index(size_t numLinks, size_t size, Varray<size_t> const &indices)
{
  std::vector<size_t> isum(size, 0);

  for (size_t n = 0; n < numLinks; ++n) { isum[indices[n]]++; }

  size_t maxIndex = 0;
  for (size_t i = 0; i < size; ++i)
    if (isum[i] > maxIndex) { maxIndex = isum[i]; }

  return maxIndex;
}

static size_t
binary_search_int(Varray<size_t> const &array, size_t len, size_t value)
{
  int64_t low = 0, high = len - 1;

  while (low <= high)
  {
    auto midpoint = low + (high - low) / 2;
    // check to see if value is equal to item in array
    if (value == array[midpoint]) return midpoint;

    (value < array[midpoint]) ? high = midpoint - 1 : low = midpoint + 1;
  }

  // item was not found
  return len;
}

static std::pair<size_t, size_t>
get_minmax_index(size_t i, size_t numLinks, Varray<size_t> const &tgtIndices)
{
  size_t minIndex = 1, maxIndex = 0;

  auto n = binary_search_int(tgtIndices, numLinks, i);
  if (n < numLinks)
  {
    minIndex = n;

    for (n = minIndex + 1; n < numLinks; ++n)
      if (i != tgtIndices[n]) break;

    maxIndex = n;

    for (n = minIndex; n > 0; --n)
      if (i != tgtIndices[n - 1]) break;

    minIndex = n;
  }

  return std::make_pair(minIndex, maxIndex);
}

template <typename T1, typename T2>
static void
remap_laf(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double tgtMissval, size_t tgtSize, RemapVars const &rv)
{
  T2 missval = tgtMissval;
  /*
    Input:
      srcArray : array with source field to be remapped

    Output:
      tgtArray : array for remapped field on target grid
  */
  auto numLinks = rv.numLinks;
  auto numWeights = rv.numWeights;             // num of weights used in remapping
  auto const &weights = rv.weights;            // remapping weights for each link
  auto const &tgtIndices = rv.tgtCellIndices;  // target indices for each link
  auto const &srcIndices = rv.srcCellIndices;  // source indices for each link

  assert(std::ranges::is_sorted(tgtIndices));

  std::ranges::fill(tgtArray, missval);

  if (numLinks == 0) return;

  auto max_cls = get_max_index(numLinks, tgtSize, tgtIndices);

#ifdef _OPENMP
  Varray2D<T1> src_cls2(Threading::ompNumMaxThreads, Varray<T1>(max_cls));
  Varray2D<double> src_weights2(Threading::ompNumMaxThreads, Varray<double>(max_cls));
#else
  Varray<T1> src_cls(max_cls);
  Varray<double> src_weights(max_cls);
#endif

  for (size_t n = 0; n < numLinks; ++n)
    if (fp_is_equal(tgtArray[tgtIndices[n]], missval)) { tgtArray[tgtIndices[n]] = 0.0; }

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic, 1)
#endif
  for (size_t i = 0; i < tgtSize; ++i)
  {
#ifdef _OPENMP
    auto ompthID = cdo_omp_get_thread_num();
    auto &src_cls = src_cls2[ompthID];
    auto &src_weights = src_weights2[ompthID];
#endif
    std::ranges::fill(src_cls, static_cast<T1>(0.0));
    std::ranges::fill(src_weights, 0.0);

    auto [minIndex, maxIndex] = get_minmax_index(i, numLinks, tgtIndices);

    size_t ncls = 0;
    for (auto j = minIndex; j < maxIndex; ++j)
    {
      auto value = srcArray[srcIndices[j]];

      size_t k = 0;
      for (; k < ncls; ++k)
        if (is_equal(value, src_cls[k])) break;

      if (k == ncls)
      {
        src_cls[k] = value;
        ncls++;
      }

      src_weights[k] += weights[numWeights * j];
    }
    // printf("i, minIndex, maxIndex, ncls %zu %zu %zu %zu\n", i, minIndex, maxIndex, ncls);

    if (ncls)
    {
      size_t imax = 0;
      auto weight = src_weights[0];
      for (size_t k = 1; k < ncls; ++k)
      {
        if (src_weights[k] > weight)
        {
          weight = src_weights[k];
          imax = k;
        }
      }

      // for (k = 0; k < ncls; ++k) printf(" i  k, src_weights[k],  src_cls[k] %zu %zu %g %g\n", i, k, src_weights[k],
      // src_cls[k]); printf("imax, src_weights[imax],  src_cls[imax] %zu %zu %g %g\n", i , imax, src_weights[imax],
      // src_cls[imax]);
      tgtArray[i] = src_cls[imax];
    }
  }
}

void
remap_laf(Field &field2, double missval, size_t gridsize2, RemapVars const &rv, Field const &field1)
{
  auto func = [&](auto const &v1, auto &v2) { remap_laf(v1, v2, missval, gridsize2, rv); };
  field_operation2(func, field1, field2);
}

template <typename T1, typename T2>
static void
remap_avg(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double tgtMissval, size_t tgtSize, RemapVars const &rv)
{
  T2 missval = tgtMissval;
  auto numLinks = rv.numLinks;
  auto numWeights = rv.numWeights;             // num of weights used in remapping
  auto const &weights = rv.weights;            // remapping weights for each link
  auto const &tgtIndices = rv.tgtCellIndices;  // target indices for each link
  auto const &srcIndices = rv.srcCellIndices;  // source indices for each link

  assert(std::ranges::is_sorted(tgtIndices));

  /*
  for (size_t n = 0; n < tgtSize; ++n) tgtArray[n] = missval;

  std::vector<int> count(tgtSize, 0);

#ifdef SX
#pragma cdir nodep
#endif
  for (size_t n = 0; n < numLinks; ++n)
    if (fp_is_equal(tgtArray[tgtIndices[n]], missval)) tgtArray[tgtIndices[n]] = 0.0;

  for (size_t n = 0; n < numLinks; ++n)
    {
      // printf("%5d %5d %5d %g # tgtIndices srcIndices n\n", tgtIndices[n], srcIndices[n], n, weights[numWeights*n]);
      // tgtArray[tgtIndices[n]] += srcArray[srcIndices[n]]*weights[numWeights*n];
      tgtArray[tgtIndices[n]] += srcArray[srcIndices[n]];
      count[tgtIndices[n]] += 1;
      if (src_cell_frac[srcIndices[n]] < 1.0)
        printf("%zu %zu %zu %g %g %g %g\n", n, tgtIndices[n], srcIndices[n], srcArray[srcIndices[n]], weights[numWeights * n],
               tgtArray[tgtIndices[n]], src_cell_frac[srcIndices[n]]);
    }

  for (size_t i = 0; i < tgtSize; ++i)
    {
      if (count[i] > 0) tgtArray[i] /= count[i];
    }
  */
  /*
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(dynamic, 1)
#endif
  */
  // size_t max_values = 0;
  for (size_t i = 0; i < tgtSize; ++i)
  {
    size_t nvalues = 0;
    double sum = 0.0;
    auto [minIndex, maxIndex] = get_minmax_index(i, numLinks, tgtIndices);

    // auto numIndices = (maxIndex - minIndex) + 1;
    // double lim = 0.1 / numIndices;
    constexpr double lim = 0.0;
    for (auto j = minIndex; j < maxIndex; ++j)
    {
      auto value = srcArray[srcIndices[j]];
      if (weights[numWeights * j] > lim && fp_is_not_equal(value, missval))
      {
        sum += value;
        // weights += weights[numWeights * n];
        nvalues++;
      }
    }

    tgtArray[i] = (nvalues > 0) ? sum / nvalues : missval;
    // printf("%zu %zu %g %g\n", i+1, nvalues, weights, sum);
    // max_values += nvalues;
  }
  // printf("max_values = %zu  numLinks = %zu\n", max_values, numLinks);
}

void
remap_avg(Field &field2, double missval, size_t gridsize2, RemapVars const &rv, Field const &field1)
{
  auto func = [&](auto const &v1, auto &v2) { remap_avg(v1, v2, missval, gridsize2, rv); };
  field_operation2(func, field1, field2);
}

void
remap_vars_init(RemapMethod mapType, int remapOrder, RemapVars &rv)
{
  // Determine the number of weights
  rv.numWeights = (mapType == RemapMethod::BICUBIC) ? 4 : 1;
  if (mapType == RemapMethod::CONSERV && remapOrder == 2) rv.numWeights = 3;
}

void
remap_vars_free(RemapVars &rv)
{
  if (rv.linksOffset.size() > 0) varray_free(rv.linksOffset);
  if (rv.linksPerValue.size() > 0) varray_free(rv.linksPerValue);
  varray_free(rv.srcCellIndices);
  varray_free(rv.tgtCellIndices);
  varray_free(rv.weights);
}

void
remap_vars_check_weights(RemapVars const &rv)
{
  auto numLinks = rv.numLinks;
  auto numWeights = rv.numWeights;
  auto normOpt = rv.normOpt;
  auto const &srcIndices = rv.srcCellIndices;
  auto const &tgtIndices = rv.tgtCellIndices;
  auto const &weights = rv.weights;

  for (size_t n = 0; n < numLinks; ++n)
  {
    if (weights[n * numWeights] < -0.01)
      cdo_print("Map weight < 0! grid1idx=%zu grid2idx=%zu nlink=%zu weights=%g", srcIndices[n], tgtIndices[n], n,
                weights[n * numWeights]);

    if (normOpt != NormOpt::NONE && weights[n * numWeights] > 1.01)
      cdo_print("Map weight > 1! grid1idx=%zu grid2idx=%zu nlink=%zu weights=%g", srcIndices[n], tgtIndices[n], n,
                weights[n * numWeights]);
  }
}
