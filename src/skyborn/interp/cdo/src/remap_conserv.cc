/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <atomic>

#include "cpp_lib.h"
#include "cdo_timer.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "cdo_omp.h"
#include "remap.h"
#include "remap_method_conserv.h"
#include "remap_store_link.h"
#include "progress.h"
#include <mpim_grid.h>

// 1st order conservative interpolation

static void
set_cell_coordinates_reg2d(size_t nx, size_t cellIndex, RemapGrid const *remapGrid, yac_grid_cell const &yacGridCell,
                           double *restrict x, double *restrict y)
{
  auto storeXY = (x && y);
  auto xyz = yacGridCell.coordinates_xyz;
  auto iy = cellIndex / nx;
  auto ix = cellIndex - iy * nx;
  auto const *cornerLonsReg2d = &remapGrid->cornerLonsReg2d[ix];
  auto const *cornerLatsReg2d = &remapGrid->cornerLatsReg2d[iy];
  constexpr int xi[4] = { 0, 1, 1, 0 };
  constexpr int yi[4] = { 0, 0, 1, 1 };
  for (int k = 0; k < 4; ++k)
  {
    auto lon = cornerLonsReg2d[xi[k]];
    auto lat = cornerLatsReg2d[yi[k]];
    gcLLtoXYZ(lon, lat, xyz[k]);
    if (storeXY)
    {
      x[k] = lon;
      y[k] = lat;
    }
  }
}

static void
set_cell_coordinates_unstruct(size_t cellIndex, size_t numCorners, RemapGrid const *remapGrid, yac_grid_cell const &yacGridCell,
                              double *restrict x, double *restrict y)
{
  auto storeXY = (x && y);
  auto xyz = yacGridCell.coordinates_xyz;

  auto const *cornerLons = &remapGrid->cornerLons[cellIndex * numCorners];
  auto const *cornerLats = &remapGrid->cornerLats[cellIndex * numCorners];
  for (size_t i = 0; i < numCorners; ++i) { gcLLtoXYZ(cornerLons[i], cornerLats[i], xyz[i]); }

  if (storeXY)
  {
    for (size_t k = 0; k < numCorners; ++k) { x[k] = cornerLons[k]; }
    for (size_t k = 0; k < numCorners; ++k) { y[k] = cornerLats[k]; }
  }
}

static void
set_cell_coordinates_yac(RemapGridType remapGridType, size_t cellIndex, size_t numCorners, RemapGrid const *remapGrid,
                         yac_grid_cell const &yacGridCell, double *x, double *y)
{
  auto isReg2D = (remapGridType == RemapGridType::Reg2D);
  isReg2D ? set_cell_coordinates_reg2d(remapGrid->dims[0], cellIndex, remapGrid, yacGridCell, x, y)
          : set_cell_coordinates_unstruct(cellIndex, numCorners, remapGrid, yacGridCell, x, y);
}

static void
set_cell_coordinates(RemapGridType remapGridType, size_t cellIndex, size_t numCorners, RemapGrid const *remapGrid,
                     GridCell const &gridCell)
{
  set_cell_coordinates_yac(remapGridType, cellIndex, numCorners, remapGrid, gridCell.yacGridCell, gridCell.coordinatesX,
                           gridCell.coordinatesY);
}

static void
set_coordinates_yac(size_t numCells, RemapGridType remapGridType, Varray<size_t> const &cellIndices, size_t numCorners,
                    RemapGrid const *remapGrid, Varray<yac_grid_cell> const &yacGridCell)
{
  for (size_t i = 0; i < numCells; ++i)
    set_cell_coordinates_yac(remapGridType, cellIndices[i], numCorners, remapGrid, yacGridCell[i], nullptr, nullptr);
}

static void
set_coordinates_yac(double const (*xyzCoords)[3], size_t numCells, Varray<size_t> const &cellIndices, size_t numCorners,
                    Varray<yac_grid_cell> const &yacGridCell)
{
  for (size_t i = 0; i < numCells; ++i)
  {
    auto offset = cellIndices[i] * numCorners;
    auto xyz = yacGridCell[i].coordinates_xyz;
    for (size_t k = 0; k < numCorners; ++k)
      for (size_t l = 0; l < 3; ++l) { xyz[k][l] = xyzCoords[offset + k][l]; }
  }
}

static void
gridcell_init_yac(GridCell &gridCell, size_t numCorners, enum yac_edge_type *edgeType)
{
  gridCell.yacGridCell.array_size = numCorners;
  gridCell.yacGridCell.num_corners = numCorners;
  gridCell.yacGridCell.edge_type = edgeType;
  gridCell.yacGridCell.coordinates_xyz = new double[numCorners][3];
  gridCell.coordinatesX = new double[numCorners];
  gridCell.coordinatesY = new double[numCorners];
}

static void
gridcell_free_yac(const GridCell &gridCell)
{
  delete[] gridCell.yacGridCell.coordinates_xyz;
  delete[] gridCell.coordinatesX;
  delete[] gridCell.coordinatesY;
}

static int
get_lonlat_circle_index(RemapGridType remapGridType, size_t gridsize, size_t numCorners, Varray<double> const &clon,
                        Varray<double> const &clat)
{
  int lonlatCircleIndex = -1;

  if (numCorners == 4)
  {
    if (remapGridType == RemapGridType::Reg2D) { lonlatCircleIndex = 1; }
    else
    {
      size_t incr = (gridsize < 100) ? 1 : gridsize / 30 - 1;
      size_t num_i = 0, num_eq0 = 0, num_eq1 = 0;

      for (size_t i = 0; i < gridsize; i += incr)
      {
        auto i4 = i * 4;
        num_i++;
        // clang-format off
        if (is_equal(clon[i4 + 1], clon[i4 + 2]) && is_equal(clon[i4 + 3], clon[i4 + 0]) &&
            is_equal(clat[i4 + 0], clat[i4 + 1]) && is_equal(clat[i4 + 2], clat[i4 + 3]))
        {
          num_eq1++;
        }
        else if (is_equal(clon[i4 + 0], clon[i4 + 1]) && is_equal(clon[i4 + 2], clon[i4 + 3]) &&
                 is_equal(clat[i4 + 1], clat[i4 + 2]) && is_equal(clat[i4 + 3], clat[i4 + 0]))
        {
          num_eq0++;
        }
        // clang-format on
      }

      if (num_i == num_eq1) { lonlatCircleIndex = 1; }
      if (num_i == num_eq0) { lonlatCircleIndex = 0; }
    }
  }

  // printf("lonlatCircleIndex %d\n", lonlatCircleIndex);

  return lonlatCircleIndex;
}

static int
get_lonlat_circle_index(RemapGrid const *remapGrid)
{
  auto gridsize = remapGrid->size;
  auto numCorners = remapGrid->numCorners;
  auto const &clon = remapGrid->cornerLons;
  auto const &clat = remapGrid->cornerLats;
  return get_lonlat_circle_index(remapGrid->type, gridsize, numCorners, clon, clat);
}

static void
normalize_weights(RemapGrid const *tgtGrid, RemapVars &rv)
{
  // Include centroids in weights and normalize using target cell area if requested
  auto numLinks = rv.numLinks;
  auto numWeights = rv.numWeights;

  if (rv.normOpt == NormOpt::DESTAREA)
  {
    auto const &cellArea = tgtGrid->cellArea;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
    for (size_t i = 0; i < numLinks; ++i)
    {
      auto index = rv.tgtCellIndices[i];  // current linear index for target grid cell
      auto normFactor = is_not_equal(cellArea[index], 0.0) ? 1.0 / cellArea[index] : 0.0;
      rv.weights[i * numWeights] *= normFactor;
    }
  }
  else if (rv.normOpt == NormOpt::FRACAREA)
  {
    auto const &cellFrac = tgtGrid->cellFrac;
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
    for (size_t i = 0; i < numLinks; ++i)
    {
      auto index = rv.tgtCellIndices[i];  // current linear index for target grid cell
      auto normFactor = is_not_equal(cellFrac[index], 0.0) ? 1.0 / cellFrac[index] : 0.0;
      rv.weights[i * numWeights] *= normFactor;
    }
  }
  else if (rv.normOpt == NormOpt::NONE) {}
}

static void
normalize_weights(NormOpt normOpt, double cellArea, double cellFrac, size_t numWeights, Varray<double> &weights)
{
  if (normOpt == NormOpt::DESTAREA)
  {
    auto normFactor = is_not_equal(cellArea, 0.0) ? 1.0 / cellArea : 0.0;
    for (size_t i = 0; i < numWeights; ++i) { weights[i] *= normFactor; }
  }
  else if (normOpt == NormOpt::FRACAREA)
  {
    auto normFactor = is_not_equal(cellFrac, 0.0) ? 1.0 / cellFrac : 0.0;
    for (size_t i = 0; i < numWeights; ++i) { weights[i] *= normFactor; }
  }
  else if (normOpt == NormOpt::NONE) {}
}

static void
correct_weights(double cellArea, size_t numWeights, Varray<double> &weights)
{
  for (size_t i = 0; i < numWeights; ++i) { weights[i] /= cellArea; }
  yac_correct_weights(numWeights, weights.data());
  for (size_t i = 0; i < numWeights; ++i) { weights[i] *= cellArea; }
}

#ifdef HAVE_LIB_RANGES_ZIP
#include <ranges>
static void
sort_weights_by_index_zip(size_t numWeights, Varray<size_t> &indices, Varray<double> &weights)
{
  /*
  static bool doPrint = true;
  if (doPrint)
    {
      doPrint = false;
      printf("using sort_weights_by_index_zip()\n");
    }
  */
  auto r = std::views::zip(indices, weights);
  std::sort(r.begin(), r.begin() + numWeights,
            [](auto a, auto b)
            {
              auto [ai, aw] = a;
              auto [bi, bw] = b;
              return (ai < bi);
            });
}
#else
static void
sort_weights_by_index(size_t numWeights, Varray<size_t> &indices, Varray<double> &weights)
{
  struct IndexWeightX
  {
    size_t index;
    double weight;
  };

  Varray<IndexWeightX> indexWeights(numWeights);

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
#endif

static void
sort_weights(size_t numWeights, Varray<size_t> &indices, Varray<double> &weights)
{
  if (numWeights <= 1) return;

  if (is_sorted_list(numWeights, indices.data())) return;

#ifdef HAVE_LIB_RANGES_ZIP
  sort_weights_by_index_zip(numWeights, indices, weights);
#else
  sort_weights_by_index(numWeights, indices, weights);
#endif
}
/*
static void
reg2d_bound_box(RemapGrid *remapGrid, double *gridBoundBox)
{
  auto nx = remapGrid->dims[0];
  auto ny = remapGrid->dims[1];
  auto const &reg2d_corner_lon = remapGrid->reg2d_corner_lon;
  auto const &reg2d_corner_lat = remapGrid->reg2d_corner_lat;

  gridBoundBox[0] = reg2d_corner_lat[0];
  gridBoundBox[1] = reg2d_corner_lat[ny];
  if (gridBoundBox[0] > gridBoundBox[1])
    {
      gridBoundBox[0] = reg2d_corner_lat[ny];
      gridBoundBox[1] = reg2d_corner_lat[0];
    }
  gridBoundBox[2] = reg2d_corner_lon[0];
  gridBoundBox[3] = reg2d_corner_lon[nx];
}
*/

static void
scale_cellfrac(size_t numCells, Varray<double> &cellFrac, Varray<double> const &cellArea)
{
  for (size_t i = 0; i < numCells; ++i)
    if (is_not_equal(cellArea[i], 0)) { cellFrac[i] /= cellArea[i]; }
}

static void
vec_index_weights(Varray<double> &vec, size_t numWeights, Varray<size_t> const &indices, Varray<double> const &weights)
{
  for (size_t i = 0; i < numWeights; ++i)
  {
    auto index = indices[i];
    auto weight = weights[i];
#ifndef __PGI
#ifdef _OPENMP
#pragma omp atomic
#endif
    vec[index] += weight;
#endif
  }
}

static size_t
remove_invalid_areas(size_t numSearchCells, Varray<size_t> &indices, Varray<double> &areas)
{
  size_t n = 0;
  for (size_t i = 0; i < numSearchCells; ++i)
  {
    if (areas[i] > 0.0)
    {
      indices[n] = indices[i];
      areas[n] = areas[i];
      n++;
    }
  }

  return n;
}

static size_t
remove_invalid_weights(size_t gridSize, size_t numWeights, Varray<size_t> &indices, Varray<double> &weights)
{
  size_t n = 0;
  for (size_t i = 0; i < numWeights; ++i)
  {
    auto index = (weights[i] > 0.0) ? indices[i] : gridSize;
    if (index != gridSize)
    {
      indices[n] = index;
      weights[n] = weights[i];
      n++;
    }
  }

  return n;
}

static size_t
remove_unmask_weights(Vmask const &gridMask, size_t numWeights, Varray<size_t> &indices, Varray<double> &weights)
{
  size_t n = 0;
  for (size_t i = 0; i < numWeights; ++i)
  {
    auto index = indices[i];
    /*
      Store the appropriate indices and weights.
      Also add contributions to cell areas.
      The source grid mask is the master mask.
    */
    if (gridMask[index])
    {
      indices[n] = index;
      weights[n] = weights[i];
      n++;
    }
  }

  return n;
}

static void
stat_update(size_t numSearchCells, size_t (&numSearchCellsStat)[3])
{
  numSearchCellsStat[0] += numSearchCells;
  numSearchCellsStat[1] = std::min(numSearchCellsStat[1], numSearchCells);
  numSearchCellsStat[2] = std::max(numSearchCellsStat[2], numSearchCells);
}

static size_t
remap_search_cells(RemapSearch &rsearch, bool isReg2dCell, GridCell const &gridCell, Varray<size_t> &searchIndices)
{
  if (rsearch.srcGrid->type == RemapGridType::Reg2D)
    return rsearch.gcs.reg2d->do_cellsearch(isReg2dCell, gridCell, searchIndices);
  else
    return rsearch.gcs.unstruct.do_cellsearch(isReg2dCell, gridCell, searchIndices);
}

void
remap_conserv_weights(RemapSearch &remapSearch, RemapVars &rv)
{
  auto srcGrid = remapSearch.srcGrid;
  auto tgtGrid = remapSearch.tgtGrid;

  auto doCheck = true;

  auto srcGridType = srcGrid->type;
  auto tgtGridType = tgtGrid->type;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  cdo::Progress progress;
  cdo::timer timer;

  auto srcGridSize = srcGrid->size;
  auto tgtGridSize = tgtGrid->size;

  auto srcNumCorners = srcGrid->numCorners;
  auto tgtNumCorners = tgtGrid->numCorners;

  enum yac_edge_type lonlatCircleType[]
      = { YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE };
  Varray<enum yac_edge_type> greatCircleType(std::max(srcNumCorners, tgtNumCorners), YAC_GREAT_CIRCLE_EDGE);

  auto srcEdgeType = greatCircleType.data();
  auto tgtEdgeType = greatCircleType.data();

  enum yac_cell_type tgtCellType = YAC_MIXED_CELL;

  if (srcNumCorners == 4)
  {
    auto lonlatCircleIndex = get_lonlat_circle_index(srcGrid);
    if (lonlatCircleIndex >= 0) { srcEdgeType = &lonlatCircleType[lonlatCircleIndex]; }
  }

  if (tgtNumCorners == 4)
  {
    auto lonlatCircleIndex = get_lonlat_circle_index(tgtGrid);
    if (lonlatCircleIndex >= 0)
    {
      tgtCellType = YAC_LON_LAT_CELL;
      tgtEdgeType = &lonlatCircleType[lonlatCircleIndex];
    }
  }

  Varray<GridCell> tgtGridCell2(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) gridcell_init_yac(tgtGridCell2[i], tgtNumCorners, tgtEdgeType);

  Varray<CellSearch> cellSearch2(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
  {
    cellSearch2[i].numCorners = srcNumCorners;
    cellSearch2[i].edgeType = srcEdgeType;
  }

  auto numCorners = srcNumCorners;  // num of corners of search cells

  // double srcGridBoundBox[4];
  // if (srcGridType == RemapGridType::Reg2D) reg2d_bound_box(srcGrid, srcGridBoundBox);

  std::vector<WeightLinks> weightLinks(tgtGridSize);

  std::atomic<long> numLinksPerValue{ -1 };
  std::atomic<size_t> atomicCount{ 0 };

  size_t numSearchCellsStat[3] = { 0, 100000, 0 };

  Varray<Varray<size_t>> indices2(Threading::ompNumMaxThreads);

  // Loop over target grid cells

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    auto &cellSearch = cellSearch2[ompthID];
    auto &indices = indices2[ompthID];
    auto &tgtGridCell = tgtGridCell2[ompthID];

    if (ompthID == 0 && tgtGridSize > progressMinSize) progress.update((double) atomicCount / tgtGridSize);

    weightLinks[tgtCellIndex].nlinks = 0;

    if (!tgtGrid->mask[tgtCellIndex]) continue;

    set_cell_coordinates(tgtGridType, tgtCellIndex, tgtNumCorners, tgtGrid, tgtGridCell);

    // Get search cells
    // numSearchCells = remap_search_cells(remapSearch, (tgtCellType == YAC_LON_LAT_CELL), tgtGridCell, indices);
    auto numSearchCells = remap_search_cells(remapSearch, (tgtGridType == RemapGridType::Reg2D), tgtGridCell, indices);

    if (1 && Options::cdoVerbose) { stat_update(numSearchCells, numSearchCellsStat); }

    if (0 && Options::cdoVerbose) cdo_print("tgtCellIndex %zu  numSearchCells %zu", tgtCellIndex, numSearchCells);

    if (numSearchCells == 0) continue;

    // Create search arrays

    cellSearch.realloc(numSearchCells);

    if (remapSearch.gcs.m_xyzCoords)
      set_coordinates_yac(remapSearch.gcs.m_xyzCoords, numSearchCells, indices, numCorners, cellSearch.gridCells);
    else
      set_coordinates_yac(numSearchCells, srcGridType, indices, numCorners, srcGrid, cellSearch.gridCells);

    if (tgtNumCorners < 4 || tgtCellType == YAC_LON_LAT_CELL)
      cdo_compute_overlap_info(numSearchCells, cellSearch, tgtGridCell);
    else
      cdo_compute_concave_overlap_info(numSearchCells, cellSearch, tgtGridCell);

    auto &partialWeights = cellSearch.partialAreas;

    auto numWeights = remove_invalid_areas(numSearchCells, indices, partialWeights);

    auto tgtCellArea = gridcell_area(tgtGridCell.yacGridCell);
    tgtGrid->cellArea[tgtCellIndex] = tgtCellArea;

    if (rv.normOpt == NormOpt::FRACAREA) correct_weights(tgtCellArea, numWeights, partialWeights);

    numWeights = remove_invalid_weights(srcGridSize, numWeights, indices, partialWeights);

    vec_index_weights(srcGrid->cellArea, numWeights, indices, partialWeights);

    numWeights = remove_unmask_weights(srcGrid->mask, numWeights, indices, partialWeights);

    vec_index_weights(srcGrid->cellFrac, numWeights, indices, partialWeights);

    tgtGrid->cellFrac[tgtCellIndex] = varray_sum(numWeights, partialWeights);

    store_weightlinks(1, numWeights, indices.data(), partialWeights.data(), tgtCellIndex, weightLinks);

    if (numWeights > 0)
    {
      if (numLinksPerValue == -1) { numLinksPerValue = numWeights; }
      else if (numLinksPerValue > 0 && numLinksPerValue != (long) numWeights) { numLinksPerValue = 0; }
    }
  }

  if (numLinksPerValue > 0) rv.numLinksPerValue = numLinksPerValue;

  if (Options::cdoVerbose)
  {
    cdo_print("Num search cells min,mean,max :  %zu  %3.1f  %zu", numSearchCellsStat[1],
              numSearchCellsStat[0] / (double) tgtGridSize, numSearchCellsStat[2]);
  }

  // Finished with all cells: deallocate search arrays
  for (auto ompthID = 0; ompthID < Threading::ompNumMaxThreads; ++ompthID)
  {
    cellSearch2[ompthID].free();
    gridcell_free_yac(tgtGridCell2[ompthID]);
  }

  weight_links_to_remap_links(1, tgtGridSize, weightLinks, rv);

  // Normalize weights using target cell area if requested
  normalize_weights(tgtGrid, rv);

  if (Options::cdoVerbose) cdo_print("Total number of links = %zu", rv.numLinks);

  scale_cellfrac(srcGridSize, srcGrid->cellFrac, srcGrid->cellArea);
  scale_cellfrac(tgtGridSize, tgtGrid->cellFrac, tgtGrid->cellArea);

  // Perform some error checking on final weights
  if (doCheck)
  {
    remap_check_area(srcGridSize, srcGrid->cellArea, "Source");
    remap_check_area(tgtGridSize, tgtGrid->cellArea, "Target");

    remap_vars_check_weights(rv);
  }

  if (Options::cdoVerbose) cdo_print("Cells search: %.2f seconds", timer.elapsed());
}  // remap_conserv_weights

template <typename T>
static double
conserv_remap(Varray<T> const &srcArray, size_t numWeights, Varray<double> const &weights, Varray<size_t> const &srcIndices)
{
  double tgtPoint = 0.0;
  for (size_t i = 0; i < numWeights; ++i) { tgtPoint += srcArray[srcIndices[i]] * weights[i]; }

  return tgtPoint;
}

template <typename T1, typename T2>
static void
remap_conserv(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double srcMissval, size_t numMissVals, NormOpt normOpt,
              RemapSearch &remapSearch)
{
  T1 missval = srcMissval;
  auto srcGrid = remapSearch.srcGrid;
  auto tgtGrid = remapSearch.tgtGrid;

  auto doCheck = true;

  // Variables necessary if segment manages to hit pole
  auto srcGridType = srcGrid->type;
  auto tgtGridType = tgtGrid->type;

  if (Options::cdoVerbose) cdo_print("Called %s()", __func__);

  cdo::Progress progress;
  cdo::timer timer;

  auto srcGridSize = srcGrid->size;
  auto tgtGridSize = tgtGrid->size;

  Vmask srcGridMask;
  if (numMissVals) remap_set_mask(srcArray, srcGridSize, numMissVals, srcMissval, srcGridMask);

  auto srcNumCorners = srcGrid->numCorners;
  auto tgtNumCorners = tgtGrid->numCorners;

  enum yac_edge_type lonlatCircleType[]
      = { YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE };
  Varray<enum yac_edge_type> greatCircleType(std::max(srcNumCorners, tgtNumCorners), YAC_GREAT_CIRCLE_EDGE);

  auto srcEdgeType = greatCircleType.data();
  auto tgtEdgeType = greatCircleType.data();

  enum yac_cell_type tgtCellType = YAC_MIXED_CELL;

  if (srcNumCorners == 4)
  {
    auto lonlatCircleIndex = get_lonlat_circle_index(srcGrid);
    if (lonlatCircleIndex >= 0) srcEdgeType = &lonlatCircleType[lonlatCircleIndex];
  }

  if (tgtNumCorners == 4)
  {
    auto lonlatCircleIndex = get_lonlat_circle_index(tgtGrid);
    if (lonlatCircleIndex >= 0)
    {
      tgtCellType = YAC_LON_LAT_CELL;
      tgtEdgeType = &lonlatCircleType[lonlatCircleIndex];
    }
  }

  Varray<GridCell> tgtGridCell2(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i) gridcell_init_yac(tgtGridCell2[i], tgtNumCorners, tgtEdgeType);

  Varray<CellSearch> cellSearch2(Threading::ompNumMaxThreads);
  for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
  {
    cellSearch2[i].numCorners = srcNumCorners;
    cellSearch2[i].edgeType = srcEdgeType;
  }

  auto numCorners = srcNumCorners;  // num of corners of search cells

  // double srcGridBoundBox[4];
  // if (srcGridType == RemapGridType::Reg2D) reg2d_bound_box(srcGrid, srcGridBoundBox);

  std::atomic<size_t> atomicCount{ 0 };

  size_t numSearchCellsStat[3] = { 0, 100000, 0 };

  Varray<Varray<size_t>> indices2(Threading::ompNumMaxThreads);

  // Loop over target grid cells

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t tgtCellIndex = 0; tgtCellIndex < tgtGridSize; ++tgtCellIndex)
  {
    atomicCount++;
    auto ompthID = cdo_omp_get_thread_num();
    auto &cellSearch = cellSearch2[ompthID];
    auto &indices = indices2[ompthID];
    auto &tgtGridCell = tgtGridCell2[ompthID];

    if (ompthID == 0 && tgtGridSize > progressMinSize) progress.update((double) atomicCount / tgtGridSize);

    tgtArray[tgtCellIndex] = missval;

    if (!tgtGrid->mask[tgtCellIndex]) continue;

    set_cell_coordinates(tgtGridType, tgtCellIndex, tgtNumCorners, tgtGrid, tgtGridCell);

    // Get search cells
    // numSearchCells = remap_search_cells(remapSearch, (tgtCellType == YAC_LON_LAT_CELL), tgtGridCell, indices);
    auto numSearchCells = remap_search_cells(remapSearch, (tgtGridType == RemapGridType::Reg2D), tgtGridCell, indices);

    if (1 && Options::cdoVerbose) { stat_update(numSearchCells, numSearchCellsStat); }

    if (0 && Options::cdoVerbose) cdo_print("tgtCellIndex %zu  numSearchCells %zu", tgtCellIndex, numSearchCells);

    if (numSearchCells == 0) continue;

    // Create search arrays

    cellSearch.realloc(numSearchCells);

    if (remapSearch.gcs.m_xyzCoords)
      set_coordinates_yac(remapSearch.gcs.m_xyzCoords, numSearchCells, indices, numCorners, cellSearch.gridCells);
    else
      set_coordinates_yac(numSearchCells, srcGridType, indices, numCorners, srcGrid, cellSearch.gridCells);

    if (tgtNumCorners < 4 || tgtCellType == YAC_LON_LAT_CELL)
      cdo_compute_overlap_info(numSearchCells, cellSearch, tgtGridCell);
    else
      cdo_compute_concave_overlap_info(numSearchCells, cellSearch, tgtGridCell);

    auto &partialWeights = cellSearch.partialAreas;

    auto numWeights = remove_invalid_areas(numSearchCells, indices, partialWeights);

    auto tgtCellArea = gridcell_area(tgtGridCell.yacGridCell);
    tgtGrid->cellArea[tgtCellIndex] = tgtCellArea;

    if (normOpt == NormOpt::FRACAREA) correct_weights(tgtCellArea, numWeights, partialWeights);

    numWeights = remove_invalid_weights(srcGridSize, numWeights, indices, partialWeights);
    if (srcGridMask.size() > 0) numWeights = remove_unmask_weights(srcGridMask, numWeights, indices, partialWeights);

    tgtGrid->cellFrac[tgtCellIndex] = varray_sum(numWeights, partialWeights);

    if (numWeights)
    {
      sort_weights(numWeights, indices, partialWeights);
      // Normalize weights using cell target area if requested
      normalize_weights(normOpt, tgtCellArea, tgtGrid->cellFrac[tgtCellIndex], numWeights, partialWeights);
      tgtArray[tgtCellIndex] = conserv_remap(srcArray, numWeights, partialWeights, indices);
    }
  }

  if (Options::cdoVerbose)
  {
    cdo_print("Num search cells min,mean,max :  %zu  %3.1f  %zu", numSearchCellsStat[1],
              numSearchCellsStat[0] / (double) tgtGridSize, numSearchCellsStat[2]);
  }

  // Finished with all cells: deallocate search arrays

  for (auto ompthID = 0; ompthID < Threading::ompNumMaxThreads; ++ompthID)
  {
    cellSearch2[ompthID].free();
    gridcell_free_yac(tgtGridCell2[ompthID]);
  }

  scale_cellfrac(tgtGridSize, tgtGrid->cellFrac, tgtGrid->cellArea);

  // Perform some error checking on final weights
  if (doCheck) remap_check_area(tgtGridSize, tgtGrid->cellArea, "Target");

  if (Options::cdoVerbose) cdo_print("Cells search: %.2f seconds", timer.elapsed());
}  // remap_conserv

void
remap_conserv(NormOpt normOpt, RemapSearch &remapSearch, Field const &field1, Field &field2)
{
  auto func = [&](auto const &v1, auto &v2) { remap_conserv(v1, v2, field1.missval, field1.numMissVals, normOpt, remapSearch); };
  field_operation2(func, field1, field2);
}

template <typename T>
static size_t
remove_missing_weights(Varray<T> const &srcArray, T missval, size_t numWeights, Varray<double> &partialWeights,
                       Varray<size_t> &indices)
{
  size_t n = 0;
  for (size_t i = 0; i < numWeights; ++i)
  {
    auto cellIndex = indices[i];
    if (fp_is_not_equal(srcArray[cellIndex], missval))
    {
      partialWeights[n] = partialWeights[i];
      indices[n] = cellIndex;
      n++;
    }
  }

  return n;
}

static double
sphere_segment_area(double latInRadian)
{
  return 2.0 * std::numbers::pi * (1.0 - std::cos(std::numbers::pi * 0.5 - latInRadian));
}

static double
latitude_area(double latMin, double latMax)
{
  return sphere_segment_area(latMin) - sphere_segment_area(latMax);
}

static void
calc_remap_indices(size_t gridsize1, size_t nv1, Varray<double> const &ybounds1, size_t ysize2, Varray<double> const &ybounds2,
                   Varray2D<size_t> &remapIndices)
{
  constexpr double scaleFactor = 1000000000.0;
  Varray<int> ymin1(gridsize1), ymax1(gridsize1);

#ifdef _OPENMP
#pragma omp parallel for if (gridsize1 > cdoMinLoopSize) schedule(static) default(shared)
#endif
  for (size_t i = 0; i < gridsize1; ++i)
  {
    auto minval = ybounds1[i * nv1];
    auto maxval = ybounds1[i * nv1];
    for (size_t k = 1; k < nv1; ++k)
    {
      auto val = ybounds1[i * nv1 + k];
      maxval = std::max(maxval, val);
      minval = std::min(minval, val);
    }
    ymin1[i] = scaleFactor * minval;
    ymax1[i] = scaleFactor * maxval;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i2 = 0; i2 < ysize2; ++i2)
  {
    int latBounds[2] = { static_cast<int>(scaleFactor * ybounds2[2 * i2]), static_cast<int>(scaleFactor * ybounds2[2 * i2 + 1]) };
    if (latBounds[0] > latBounds[1]) std::swap(latBounds[0], latBounds[1]);

    size_t numSearchCells = 0;
    for (size_t i1 = 0; i1 < gridsize1; ++i1)
    {
      if (ymin1[i1] < latBounds[1] && ymax1[i1] > latBounds[0]) { numSearchCells++; }
    }

    remapIndices[i2].resize(numSearchCells);
    size_t n = 0;
    for (size_t i1 = 0; i1 < gridsize1; ++i1)
    {
      if (ymin1[i1] < latBounds[1] && ymax1[i1] > latBounds[0]) { remapIndices[i2][n++] = i1; }
    }

    // printf("lat %zu found %zu of %zu\n", i2 + 1, numSearchCells, gridsize1);
  }
}

void
remap_weights_zonal_mean(int gridID1, int gridID2, Varray2D<size_t> &remapIndices, Varray2D<double> &remapWeights)
{
  auto gridsize1 = gridInqSize(gridID1);
  size_t nv1 = gridInqNvertex(gridID1);
  Varray<double> xbounds1(gridsize1 * nv1), ybounds1(gridsize1 * nv1);
  gridInqXbounds(gridID1, xbounds1.data());
  gridInqYbounds(gridID1, ybounds1.data());

  // Convert lonlat units if required
  cdo_grid_to_radian(gridID1, CDI_XAXIS, xbounds1, "source grid longitude bounds");
  cdo_grid_to_radian(gridID1, CDI_YAXIS, ybounds1, "source grid latitude bounds");

  auto ysize2 = gridInqYsize(gridID2);
  Varray<double> ybounds2(ysize2 * 2);
  gridInqYbounds(gridID2, ybounds2.data());

  // Convert lat units if required
  cdo_grid_to_radian(gridID2, CDI_YAXIS, ybounds2, "target grid latitude bounds");

  remapIndices.resize(ysize2);
  remapWeights.resize(ysize2);

  calc_remap_indices(gridsize1, nv1, ybounds1, ysize2, ybounds2, remapIndices);

  enum yac_edge_type lonlatCircleType[]
      = { YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE };
  Varray<enum yac_edge_type> greatCircleType(nv1, YAC_GREAT_CIRCLE_EDGE);

  auto srcEdgeType = greatCircleType.data();

  if (nv1 == 4)
  {
    auto lonlatCircleIndex = get_lonlat_circle_index(RemapGridType::Undefined, gridsize1, nv1, xbounds1, ybounds1);
    if (lonlatCircleIndex >= 0) srcEdgeType = &lonlatCircleType[lonlatCircleIndex];
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i2 = 0; i2 < ysize2; ++i2)
  {
    auto normOpt(NormOpt::FRACAREA);
    double latBounds[2] = { ybounds2[2 * i2], ybounds2[2 * i2 + 1] };
    if (latBounds[0] > latBounds[1]) std::swap(latBounds[0], latBounds[1]);

    auto tgtCellArea = latitude_area(latBounds[0], latBounds[1]);
    // printf("tgtCellArea %zu %g\n", i2 + 1, tgtCellArea);

    auto numSearchCells = remapIndices[i2].size();

    CellSearch cellSearch;
    cellSearch.numCorners = nv1;
    cellSearch.edgeType = srcEdgeType;

    cellSearch.realloc(numSearchCells);

    for (size_t j = 0; j < numSearchCells; ++j)
    {
      auto cellIndex = remapIndices[i2][j];
      auto xyz = cellSearch.gridCells[j].coordinates_xyz;
      auto const *cell_corner_lon = &xbounds1[cellIndex * nv1];
      auto const *cell_corner_lat = &ybounds1[cellIndex * nv1];
      for (size_t i = 0; i < nv1; ++i) gcLLtoXYZ(cell_corner_lon[i], cell_corner_lat[i], xyz[i]);
    }

    auto &partialWeights = cellSearch.partialAreas;
    auto &overlapCells = cellSearch.overlapCells;

    // Do the clipping and get the cell for the overlapping area
    yac_cell_lat_clipping(numSearchCells, cellSearch.gridCells.data(), latBounds, overlapCells.data());

    // Get the partial areas for the overlapping regions
    for (size_t i = 0; i < numSearchCells; ++i) { partialWeights[i] = gridcell_area(overlapCells[i]); }

    auto numWeights = remove_invalid_areas(numSearchCells, remapIndices[i2], partialWeights);
    // printf("numWeights: %zu %zu\n", numSearchCells, numWeights);

    if (normOpt == NormOpt::FRACAREA) correct_weights(tgtCellArea, numWeights, partialWeights);

    numWeights = remove_invalid_weights(gridsize1, numWeights, remapIndices[i2], partialWeights);

    remapWeights[i2].resize(numWeights);
    for (size_t i = 0; i < numWeights; ++i) remapWeights[i2][i] = partialWeights[i];

    cellSearch.free();
  }
}

template <typename T1, typename T2>
static size_t
remap_zonal_mean(Varray<T1> const &srcArray, Varray<T2> &tgtArray, double srcMissval, Varray2D<size_t> const &remapIndices,
                 Varray2D<double> const &remapWeights)
{
  T1 missval = srcMissval;
  auto ysize2 = remapIndices.size();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (size_t i2 = 0; i2 < ysize2; ++i2)
  {
    auto numWeights = remapWeights[i2].size();
    Varray<size_t> indices(numWeights);
    for (size_t i = 0; i < numWeights; ++i) { indices[i] = remapIndices[i2][i]; }
    Varray<double> partialWeights(numWeights);
    for (size_t i = 0; i < numWeights; ++i) { partialWeights[i] = remapWeights[i2][i]; }

    numWeights = remove_missing_weights(srcArray, missval, numWeights, partialWeights, indices);

    tgtArray[i2] = missval;

    if (numWeights)
    {
      auto normOpt(NormOpt::FRACAREA);
      auto tgtCellArea = 0.0;  // not needed for NormOpt::FRACAREA
      auto tgtCellFrac = varray_sum(numWeights, partialWeights);
      sort_weights(numWeights, indices, partialWeights);
      // Normalize weights using cell target area if requested
      normalize_weights(normOpt, tgtCellArea, tgtCellFrac, numWeights, partialWeights);
      tgtArray[i2] = conserv_remap(srcArray, numWeights, partialWeights, indices);
    }
  }

  size_t numMissVals = 0;
  for (size_t i2 = 0; i2 < ysize2; ++i2)
    if (fp_is_equal(tgtArray[i2], missval)) { numMissVals++; }

  return numMissVals;
}

void
remap_zonal_mean(Varray2D<size_t> const &remapIndices, Varray2D<double> const &remapWeights, Field const &field1, Field &field2)
{
  auto func = [&](auto const &v1, auto &v2) { remap_zonal_mean(v1, v2, field1.missval, remapIndices, remapWeights); };
  field_operation2(func, field1, field2);
}
