/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>

#include "cpp_lib.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_timer.h"
#include "cdo_omp.h"
#include "remap.h"
#include "remap_store_link.h"
#include "progress.h"
#include "grid_healpix.h"

// bilinear interpolation

static inline void
limit_dphi_bounds(double &dphi)
{
  if (dphi > 3.0 * PIH) { dphi -= PI2; }
  else if (dphi < -3.0 * PIH) { dphi += PI2; }
}

std::pair<double, double>
remap_find_weights(PointLonLat const &pointLL, const SquareCorners &squareCorners)
{
  constexpr double converge = 1.0e-10;  // Convergence criterion
  constexpr long MaxIteration = 100;    // Max iteration count for i, j iteration

  auto const &lons = squareCorners.lons;
  auto const &lats = squareCorners.lats;

  // Iterate to find xfrac,yfrac for bilinear approximation

  // some latitude  differences
  auto dth1 = lats[1] - lats[0];
  auto dth2 = lats[3] - lats[0];
  auto dth3 = lats[2] - lats[1] - dth2;

  // some longitude differences
  auto dph1 = lons[1] - lons[0];
  auto dph2 = lons[3] - lons[0];
  auto dph3 = lons[2] - lons[1];

  limit_dphi_bounds(dph1);
  limit_dphi_bounds(dph2);
  limit_dphi_bounds(dph3);

  dph3 -= dph2;

  // current guess for bilinear coordinate
  double xguess = 0.5;
  double yguess = 0.5;

  long iter = 0;  // iteration counters
  for (iter = 0; iter < MaxIteration; ++iter)
  {
    auto dthp = pointLL.lat() - lats[0] - dth1 * xguess - dth2 * yguess - dth3 * xguess * yguess;
    auto dphp = pointLL.lon() - lons[0];

    limit_dphi_bounds(dphp);

    dphp -= dph1 * xguess + dph2 * yguess + dph3 * xguess * yguess;

    auto mat1 = dth1 + dth3 * yguess;
    auto mat2 = dth2 + dth3 * xguess;
    auto mat3 = dph1 + dph3 * yguess;
    auto mat4 = dph2 + dph3 * xguess;

    auto determinant = mat1 * mat4 - mat2 * mat3;

    auto deli = (dthp * mat4 - dphp * mat2) / determinant;
    auto delj = (dphp * mat1 - dthp * mat3) / determinant;

    if (std::fabs(deli) < converge && std::fabs(delj) < converge) break;

    xguess += deli;
    yguess += delj;
  }

  if (iter >= MaxIteration) xguess = yguess = -1.0;

  return std::make_pair(xguess, yguess);
}

static void
bilinear_set_weights(double xfrac, double yfrac, double (&weights)[4])
{
  // clang-format off
  weights[0] = (1.0 - xfrac) * (1.0 - yfrac);
  weights[1] =        xfrac  * (1.0 - yfrac);
  weights[2] =        xfrac  *        yfrac;
  weights[3] = (1.0 - xfrac) *        yfrac;
  // clang-format on
}

int
num_src_points(Vmask const &mask, size_t const (&indices)[4], double (&lats)[4])
{
  int num = 4;

  for (int i = 0; i < 4; ++i)
  {
    if (mask[indices[i]] == 0)
    {
      num--;
      lats[i] = 0.0;
    }
  }

  return num;
}

static void
renormalize_weights(const double (&srcLats)[4], double (&weights)[4])
{
  // sum of weights for normalization
  auto sumWeights = std::fabs(srcLats[0]) + std::fabs(srcLats[1]) + std::fabs(srcLats[2]) + std::fabs(srcLats[3]);
  for (int i = 0; i < 4; ++i) weights[i] = std::fabs(srcLats[i]) / sumWeights;
}

#ifdef HAVE_LIB_RANGES_ZIP
#include <ranges>
static void
bilinear_sort_weights_by_index_zip(size_t (&indices)[4], double (&weights)[4])
{
  auto r = std::views::zip(indices, weights);
  std::sort(r.begin(), r.end(),
            [](auto a, auto b)
            {
              auto [ai, aw] = a;
              auto [bi, bw] = b;
              return (ai < bi);
            });
}
#endif

static void
bilinear_sort_weights_by_index(size_t (&indices)[4], double (&weights)[4])
{
  constexpr size_t numWeights = 4;

  struct IndexWeightX
  {
    size_t index;
    double weight;
  };

  std::array<IndexWeightX, numWeights> indexWeights;

  for (size_t i = 0; i < numWeights; ++i)
  {
    indexWeights[i].index = indices[i];
    indexWeights[i].weight = weights[i];
  }

  std::ranges::sort(indexWeights, {}, &IndexWeightX::index);

  for (size_t i = 0; i < numWeights; ++i)
  {
    indices[i] = indexWeights[i].index;
    weights[i] = indexWeights[i].weight;
  }
}

static void
bilinear_sort_weights(size_t (&indices)[4], double (&weights)[4])
{
  constexpr size_t numWeights = 4;
  if (is_sorted_list(numWeights, indices)) return;

#ifdef HAVE_LIB_RANGES_ZIP
  bilinear_sort_weights_by_index_zip(indices, weights);
#else
  bilinear_sort_weights_by_index(indices, weights);
#endif
}

static void
bilinear_warning()
{
  static auto printWarning = true;
  if (Options::cdoVerbose || printWarning)
  {
    printWarning = false;
    cdo_warning("Bilinear interpolation failed for some grid points - used a distance-weighted average instead!");
  }
}

static void
remap_bilinear_weights_regular(RemapSearch &rsearch, Vmask const &srcGridMask, PointLonLat const &pointLL, double &tgtCellFrac,
                               size_t tgtCellIndex, std::vector<WeightLinks> &weightLinks)
{
  SquareCorners squareCorners;
  double weights[4];  //  bilinear weights for four corners

  // Find nearest square of grid points on source grid
  auto searchResult = remap_search_square(rsearch, pointLL, squareCorners);

  // Check to see if points are mask points
  if (searchResult > 0) searchResult = remap_check_mask_indices(squareCorners.indices, srcGridMask);

  // If point found, find local xfrac/yfrac coordinates for weights
  if (searchResult > 0)
  {
    tgtCellFrac = 1.0;

    auto [xfrac, yfrac] = remap_find_weights(pointLL, squareCorners);
    if (xfrac >= 0.0 && yfrac >= 0.0)
    {
      // Successfully found xfrac, yfrac - compute weights
      bilinear_set_weights(xfrac, yfrac, weights);
      store_weightlinks(0, 4, squareCorners.indices, weights, tgtCellIndex, weightLinks);
    }
    else
    {
      bilinear_warning();
      searchResult = -1;
    }
  }

  // Search for bilinear failed - use a distance-weighted average instead
  // (this is typically near the pole) Distance was stored in srcLats!
  if (searchResult < 0)
  {
    if (num_src_points(srcGridMask, squareCorners.indices, squareCorners.lats) > 0)
    {
      tgtCellFrac = 1.0;
      renormalize_weights(squareCorners.lats, weights);
      store_weightlinks(0, 4, squareCorners.indices, weights, tgtCellIndex, weightLinks);
    }
  }
}

static void
remap_bilinear_weights_healpix(RemapSearch const &rsearch, Vmask const &srcGridMask, PointLonLat const &pointLL,
                               double &tgtCellFrac, size_t tgtCellIndex, std::vector<WeightLinks> &weightLinks)
{
  double weights[4];  // bilinear weights for four corners
  size_t indices[4];  // indices for the four source points

  hp_bilinear_interpolate_weights(rsearch.srcGrid->hpParams, pointLL.lon(), pointLL.lat(), indices, weights);

  // Check to see if points are mask points
  auto searchResult = remap_check_mask_indices(indices, srcGridMask);
  if (searchResult > 0)
  {
    tgtCellFrac = 1.0;
    bilinear_sort_weights(indices, weights);
    store_weightlinks(0, 4, indices, weights, tgtCellIndex, weightLinks);
  }
}

// -----------------------------------------------------------------------
// This routine computes the weights for a bilinear interpolation.
// -----------------------------------------------------------------------
void
remap_bilinear_weights(RemapSearch &rsearch, RemapVars &rv)
{
  const auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  auto isHealpixGrid = (srcGrid->type == RemapGridType::HealPix);

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (!isHealpixGrid && srcGrid->rank != 2)
    cdo_abort("Can't do bilinear interpolation if the source grid is not a regular 2D grid!");

  cdo::timer timer;
  cdo::Progress progress;

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;

  std::vector<WeightLinks> weightLinks(tgtGridSize);
  weight_links_alloc(4, tgtGridSize, weightLinks);

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && tgtGridSize > progressMinSize) { progress.update((double) atomicCount / tgtGridSize); }

    weightLinks[tgtCellIndex].nlinks = 0;

    if (!tgtGrid->mask[tgtCellIndex]) continue;

    auto pointLL = tgtGrid->get_lonlat(tgtCellIndex);

    auto &tgtCellFrac = tgtGrid->cellFrac[tgtCellIndex];
    if (isHealpixGrid)
      remap_bilinear_weights_healpix(rsearch, srcGrid->mask, pointLL, tgtCellFrac, tgtCellIndex, weightLinks);
    else
      remap_bilinear_weights_regular(rsearch, srcGrid->mask, pointLL, tgtCellFrac, tgtCellIndex, weightLinks);
  }

  weight_links_to_remap_links(0, tgtGridSize, weightLinks, rv);

  rv.numLinksPerValue = 4;

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, timer.elapsed());
}  // remap_bilinear_weights

template <typename T>
static inline T
bilinear_remap(Varray<T> const &srcArray, size_t const (&ind)[4], double const (&wgt)[4])
{
  return srcArray[ind[0]] * wgt[0] + srcArray[ind[1]] * wgt[1] + srcArray[ind[2]] * wgt[2] + srcArray[ind[3]] * wgt[3];
}

template <typename T1, typename T2>
static void
remap_bilinear_regular(RemapSearch &rsearch, Varray<T1> const &srcArray, Vmask const &srcGridMask, PointLonLat const &pointLL,
                       T2 &tgtValue)
{
  SquareCorners squareCorners;
  double weights[4];  //  bilinear weights for four corners

  // Find nearest square of grid points on source grid
  auto searchResult = remap_search_square(rsearch, pointLL, squareCorners);

  // Check to see if points are mask points
  if (searchResult > 0) searchResult = remap_check_mask_indices(squareCorners.indices, srcGridMask);

  // If point found, find local xfrac/yfrac coordinates for weights
  if (searchResult > 0)
  {
    auto [xfrac, yfrac] = remap_find_weights(pointLL, squareCorners);
    if (xfrac >= 0.0 && yfrac >= 0.0)
    {
      // Successfully found xfrac, yfrac - compute weights
      bilinear_set_weights(xfrac, yfrac, weights);
      bilinear_sort_weights(squareCorners.indices, weights);
      tgtValue = bilinear_remap(srcArray, squareCorners.indices, weights);
    }
    else
    {
      bilinear_warning();
      searchResult = -1;
    }
  }

  // Search for bilinear failed - use a distance-weighted average instead
  // (this is typically near the pole) Distance was stored in srcLats!
  if (searchResult < 0)
  {
    if (srcGridMask.size() == 0 || num_src_points(srcGridMask, squareCorners.indices, squareCorners.lats) > 0)
    {
      renormalize_weights(squareCorners.lats, weights);
      bilinear_sort_weights(squareCorners.indices, weights);
      tgtValue = bilinear_remap(srcArray, squareCorners.indices, weights);
    }
  }
}

template <typename T1, typename T2>
static void
remap_bilinear_healpix(const RemapSearch &rsearch, Varray<T1> const &srcArray, Vmask const &srcGridMask, PointLonLat const &pointLL,
                       T2 &tgtValue)
{
  size_t indices[4];  // indices for the four source points
  double weights[4];  // bilinear weights for four corners

  hp_bilinear_interpolate_weights(rsearch.srcGrid->hpParams, pointLL.lon(), pointLL.lat(), indices, weights);

  // Check to see if points are mask points
  auto searchResult = remap_check_mask_indices(indices, srcGridMask);
  if (searchResult > 0)
  {
    bilinear_sort_weights(indices, weights);
    tgtValue = bilinear_remap(srcArray, indices, weights);
  }
}

// -----------------------------------------------------------------------
// This routine computes and apply the weights for a bilinear interpolation.
// -----------------------------------------------------------------------
template <typename T1, typename T2>
static void
remap_bilinear(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double srcMissval, size_t numMissVals, RemapSearch &rsearch)
{
  T1 missval = srcMissval;
  const auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  auto isHealpixGrid = (srcGrid->type == RemapGridType::HealPix);

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  if (!isHealpixGrid && srcGrid->rank != 2)
    cdo_abort("Can't do bilinear interpolation if the source grid is not a regular 2D grid!");

  cdo::timer timer;
  cdo::Progress progress;

  auto tgtGridSize = tgtGrid->size;
  auto srcGridSize = srcGrid->size;

  Vmask srcGridMask;
  if (numMissVals) remap_set_mask(srcArray, srcGridSize, numMissVals, srcMissval, srcGridMask);

  // Compute mappings from source to target grid

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && tgtGridSize > progressMinSize) progress.update((double) atomicCount / tgtGridSize);

    auto &tgtValue = tgtArray[tgtCellIndex];
    tgtValue = missval;

    if (!tgtGrid->mask[tgtCellIndex]) continue;

    auto pointLL = tgtGrid->get_lonlat(tgtCellIndex);

    if (isHealpixGrid)
      remap_bilinear_healpix(rsearch, srcArray, srcGridMask, pointLL, tgtValue);
    else
      remap_bilinear_regular(rsearch, srcArray, srcGridMask, pointLL, tgtValue);
  }

  if (Options::cdoVerbose) cdo_print("%s: %.2f seconds", __func__, timer.elapsed());
}  // remap_bilinear

void
remap_bilinear(RemapSearch &remapSearch, Field const &field1, Field &field2)
{
  auto func = [&](auto const &v1, auto &v2) { remap_bilinear(v1, v2, field1.missval, field1.numMissVals, remapSearch); };
  field_operation2(func, field1, field2);
}
