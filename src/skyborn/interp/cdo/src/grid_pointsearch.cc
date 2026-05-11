/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstdio>

#include "pointsearch_nanoflann.h"
#include "pointsearch_spherepart.h"
#include "pointsearch_kdtree.h"
#include "pointsearch_full.h"
#include "grid_pointsearch.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include <mpim_grid.h>

UnstructMethod unstructMethod(UnstructMethod::undefined);

void
set_pointsearch_method(std::string const &methodStr)
{
  // clang-format off
  if      (methodStr == "kdtree")     unstructMethod = UnstructMethod::kdtree;
  else if (methodStr == "nanoflann")  unstructMethod = UnstructMethod::nanoflann;
  else if (methodStr == "spherepart") unstructMethod = UnstructMethod::spherepart;
  else if (methodStr == "full")       unstructMethod = UnstructMethod::full;
  else cdo_abort("Grid point search method %s not available!", methodStr);
  // clang-format on
}

static const char *
get_methodStr(UnstructMethod method)
{
  if (method == UnstructMethod::kdtree) return "kdtree";
  if (method == UnstructMethod::nanoflann) return "nanoflann";
  if (method == UnstructMethod::spherepart) return "spherepart";
  if (method == UnstructMethod::full) return "full";
  return "unknown";
}

static void
grid_pointsearch_create_healpix(GridPointsearch &gps, const RemapGrid &remapGrid)
{
  gps.healpix = new PointsearchHealpix(remapGrid.hpParams);
}

static void
grid_pointsearch_create_reg2d(GridPointsearch &gps, RemapGrid const &remapGrid)
{
  gps.reg2d = new PointsearchReg2d(remapGrid.centerLonsReg2d, remapGrid.centerLatsReg2d, gps.params);
}

static void
pointsearch_create_unstruct(GridPointsearch &gps, Varray<double> const &lons, Varray<double> const &lats)
{
  auto numPoints = lons.size();
  if (numPoints == 0) return;

  gps.numPoints = numPoints;

  gps.plons = lons.data();
  gps.plats = lats.data();

  auto &params = gps.params;
  params.isCurve = (numPoints != 1 && numPoints == params.dims[0] * params.dims[1]);

  if (unstructMethod != UnstructMethod::undefined) { gps.unstructMethod = unstructMethod; }
  if (Options::cdoVerbose) cdo_print("Point search method: %s", get_methodStr(gps.unstructMethod));

  auto method = gps.unstructMethod;
  // clang-format off
  if      (method == UnstructMethod::kdtree)     gps.unstruct.set_strategy(new PointsearchKDtree(lons, lats, params));
  else if (method == UnstructMethod::nanoflann)  gps.unstruct.set_strategy(new PointsearchNanoflann(lons, lats, params));
  else if (method == UnstructMethod::spherepart) gps.unstruct.set_strategy(new PointsearchSpherepart(lons, lats, params));
  else if (method == UnstructMethod::full)       gps.unstruct.set_strategy(new PointsearchFull(lons, lats, params));
  /*
  if      (method == UnstructMethod::kdtree)     gps.unstruct.set_strategy(std::make_unique<PointsearchKDtree>(lons, lats, params));
  else if (method == UnstructMethod::nanoflann)  gps.unstruct.set_strategy(std::make_unique<PointsearchNanoflann>(lons, lats, params));
  else if (method == UnstructMethod::spherepart) gps.unstruct.set_strategy(std::make_unique<PointsearchSpherepart>(lons, lats, params));
  else if (method == UnstructMethod::full)       gps.unstruct.set_strategy(std::make_unique<PointsearchFull>(lons, lats, params));
  */
  else cdo_abort("%s::method undefined!", __func__);
  // clang-format on
}

void
grid_pointsearch_create(GridPointsearch &gps, const RemapGrid &remapGrid)
{
  auto &params = gps.params;
  params.dims[0] = remapGrid.dims[0];
  params.dims[1] = remapGrid.dims[1];
  params.isCyclic = remapGrid.isCyclic;

  if (remapGrid.type == RemapGridType::HealPix) { grid_pointsearch_create_healpix(gps, remapGrid); }
  else if (remapGrid.type == RemapGridType::Reg2D) { grid_pointsearch_create_reg2d(gps, remapGrid); }
  else
  {
    params.useBoundBox = true;
    pointsearch_create_unstruct(gps, remapGrid.centerLons, remapGrid.centerLats);
  }
}

void
grid_pointsearch_create_unstruct(GridPointsearch &gps, Varray<double> const &lons, Varray<double> const &lats, bool useBoundBox)
{
  auto &params = gps.params;
  params.dims[0] = lons.size();
  params.dims[1] = 0;
  params.useBoundBox = useBoundBox;

  pointsearch_create_unstruct(gps, lons, lats);
}

static size_t
llindex_in_quad(size_t nx, size_t ny, size_t index, PointLonLat const &pointLL, double const *centerLons, double const *centerLats,
                bool isCyclic)
{
  if (index != GPS_NOT_FOUND)
  {
    SquareCorners squareCorners;
    for (int k = 0; k < 4; ++k)
    {
      // Determine neighbor addresses
      auto j = index / nx;
      auto i = index - j * nx;
      if (k == 1 || k == 3) i = (i > 0) ? i - 1 : (isCyclic ? nx - 1 : 0);
      if (k == 2 || k == 3) j = (j > 0) ? j - 1 : 0;

      if (point_in_quad(isCyclic, nx, ny, i, j, squareCorners, pointLL.lon(), pointLL.lat(), centerLons, centerLats)) return index;
    }
  }

  return GPS_NOT_FOUND;
}

size_t
grid_pointsearch_nearest(GridPointsearch &gps, PointLonLat const &pointLL, size_t *index, double *dist)
{
  auto numIndices = gps.unstruct.search_nearest(pointLL, index, dist);
  if (numIndices > 0)
  {
    auto const &params = gps.params;
    auto indexResult = *index;
    if (!params.extrapolation && params.isCurve)
      indexResult = llindex_in_quad(params.dims[0], params.dims[1], *index, pointLL, gps.plons, gps.plats, params.isCyclic);
    if (indexResult != GPS_NOT_FOUND) return 1;
  }

  return 0;
}

static void
select_points_in_cell(GridPointsearch &gps, PointLonLat const &pointLL, size_t *indices, double *dist, size_t &numIndices)
{
  auto maxIndices = numIndices;
  numIndices = 0;
  for (size_t i = 0; i < maxIndices; ++i)
  {
    auto const &params = gps.params;
    auto index = llindex_in_quad(params.dims[0], params.dims[1], indices[i], pointLL, gps.plons, gps.plats, params.isCyclic);
    if (index != GPS_NOT_FOUND)
    {
      indices[numIndices] = indices[i];
      dist[numIndices] = dist[i];
      numIndices++;
    }
  }
}

size_t
grid_pointsearch_qnearest(GridPointsearch &gps, PointLonLat const &pointLL, size_t nnn, size_t *indices, double *dist)
{
  auto numIndices = gps.unstruct.search_qnearest(pointLL, nnn, indices, dist);

  auto const &params = gps.params;
  if (!params.extrapolation && params.isCurve) select_points_in_cell(gps, pointLL, indices, dist, numIndices);

  return numIndices;
}
