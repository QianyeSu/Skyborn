/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>

#include "cdo_timer.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "remap.h"
#include "remap_store_link.h"
#include "progress.h"
#include "cdo_omp.h"

// -----------------------------------------------------------------------
// This routine computes the weights for a k-nearest-neighbor interpolation
// -----------------------------------------------------------------------
void
remap_knn_weights(KnnParams const &knnParams, RemapSearch &rsearch, RemapVars &rv)
{
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  cdo::Progress progress;

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;

  std::vector<WeightLinks> weightLinks(tgtGridSize);
  weight_links_alloc(knnParams.k, tgtGridSize, weightLinks);

  std::vector<KnnData> knnDataList;
  knnDataList.reserve(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) knnDataList.emplace_back(knnParams);

  cdo::timer timer;

  std::atomic<long> numLinksPerValue{ -1 };
  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && tgtGridSize > progressMinSize) progress.update((double) atomicCount / tgtGridSize);

    weightLinks[tgtCellIndex].nlinks = 0;

    if (!tgtGrid->mask[tgtCellIndex]) continue;

    auto &knnData = knnDataList[ompthID];
    auto pointLL = tgtGrid->get_lonlat(tgtCellIndex);

    // Find nearest grid points on source grid and distances to each point
    remap_search_points(rsearch, pointLL, knnData);

    // Compute weights if mask is false, eliminate those points
    auto numWeights = knnData.compute_weights(srcGrid->mask);
    if (numWeights) tgtGrid->cellFrac[tgtCellIndex] = 1.0;

    // Store the link
    store_weightlinks(0, numWeights, knnData.m_indices.data(), knnData.m_dist.data(), tgtCellIndex, weightLinks);

    if (knnParams.k > 1 && numWeights > 0)
    {
      if (numLinksPerValue == -1)
        numLinksPerValue = numWeights;
      else if (numLinksPerValue > 0 && numLinksPerValue != (long) numWeights)
        numLinksPerValue = 0;
    }
  }

  weight_links_to_remap_links(0, tgtGridSize, weightLinks, rv);

  if (knnParams.k == 1)
    rv.numLinksPerValue = 1;
  else if (numLinksPerValue > 0)
    rv.numLinksPerValue = numLinksPerValue;

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", timer.elapsed());
}  // remap_knn_weights

// -----------------------------------------------------------------------
// This routine computes and apply weights for a k-nearest-neighbor interpolation
// -----------------------------------------------------------------------
template <typename T1, typename T2>
static void
remap_knn(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double srcMissval, size_t numMissVals, KnnParams const &knnParams,
          RemapSearch &rsearch)
{
  T1 missval = srcMissval;
  auto srcGrid = rsearch.srcGrid;
  auto tgtGrid = rsearch.tgtGrid;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  cdo::Progress progress;

  // Compute mappings from source to target grid

  auto tgtGridSize = tgtGrid->size;
  auto srcGridSize = srcGrid->size;

  Vmask srcGridMask;
  if (numMissVals) remap_set_mask(srcArray, srcGridSize, numMissVals, srcMissval, srcGridMask);

  std::vector<KnnData> knnDataList;
  knnDataList.reserve(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) knnDataList.emplace_back(knnParams);

  cdo::timer timer;

  std::atomic<size_t> atomicCount{ 0 };

  // Loop over target grid
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    if (ompthID == 0 && tgtGridSize > progressMinSize) progress.update((double) atomicCount / tgtGridSize);

    auto &tgtValue = tgtArray[tgtCellIndex];
    tgtValue = missval;

    if (tgtGrid->mask.size() && !tgtGrid->mask[tgtCellIndex]) continue;

    auto &knnData = knnDataList[ompthID];
    auto pointLL = tgtGrid->get_lonlat(tgtCellIndex);

    // Find nearest grid points on source grid and distances to each point
    remap_search_points(rsearch, pointLL, knnData);

    // Compute weights if mask is false, eliminate those points
    auto numWeights = (srcGridMask.size() > 0) ? knnData.compute_weights(srcGridMask) : knnData.compute_weights();
    if (numWeights) tgtValue = knnData.array_weights_sum(srcArray);
  }

  if (Options::cdoVerbose) cdo_print("Point search nearest: %.2f seconds", timer.elapsed());
}  // remap_knn

void
remap_knn(KnnParams const &knnParams, RemapSearch &remapSearch, Field const &field1, Field &field2)
{
  auto func = [&](auto const &v1, auto &v2) { remap_knn(v1, v2, field1.missval, field1.numMissVals, knnParams, remapSearch); };
  field_operation2(func, field1, field2);
}

template <typename T1, typename T2>
void
intgrid_knn(Varray<T1> const &srcArray, Varray<T2> &tgtArray, int gridID1, int gridID2, double srcMissval, size_t numMissVals,
            KnnParams const &knnParams)
{
  auto mapType = RemapMethod::KNN;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  RemapType remap;
  remap_set_option(RemapOption::GenerateWeights, 0);
  remap_init_grids(mapType, knnParams.extrapolate, gridID1, remap.srcGrid, gridID2, remap.tgtGrid);
  remap_search_init(mapType, remap.search, remap.srcGrid, remap.tgtGrid);

  remap_knn(srcArray, tgtArray, srcMissval, numMissVals, knnParams, remap.search);

  remap_grid_free(remap.srcGrid);
  remap_grid_free(remap.tgtGrid);
}  // intgrid_knn

void
intgrid_knn(KnnParams const &knnParams, Field const &field1, Field &field2)
{
  auto func = [&](auto const &v1, auto &v2)
  { intgrid_knn(v1, v2, field1.grid, field2.grid, field1.missval, field1.numMissVals, knnParams); };
  field_operation2(func, field1, field2);

  field2.numMissVals = field_num_mv(field2);
}

void
intgrid_1nn(Field const &field1, Field &field2)
{
  KnnParams knnParams;
  knnParams.k = 1;
  intgrid_knn(knnParams, field1, field2);
}
