/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "remap.h"
#include <mpim_grid.h>
#include "cdo_output.h"

static void
grid_search_point_healpix(GridPointsearch &gps, PointLonLat const &pointLL, KnnData &knnData)
{
  /*
    Input variables:

      pointLL : longitude/latitude of the search point

    Output variables:

      knnData.m_indices[numNeighbors] : index of each of the closest points
      knnData.m_dist[numNeighbors] : distance to each of the closest points
  */
  knnData.m_numNeighbors = 0;
  auto numNeighbors = knnData.maxNeighbors();
  if (numNeighbors > 9) cdo_abort("Max number of neighbors is 9!");

  // Initialize distance and index arrays
  knnData.init_indices();
  knnData.init_dist();

  gps.healpix->search(pointLL, knnData, gps.params.searchRadius);
}

// This routine finds the closest numNeighbor points to a search point and computes a distance to each of the neighbors
static void
grid_search_point_reg2d(GridPointsearch &gps, PointLonLat const &pointLL, KnnData &knnData)
{
  /*
    Input variables:

      pointLL : longitude/latitude of the search point

    Output variables:

      knnData.m_indices[numNeighbors] : index of each of the closest points
      knnData.m_dist[numNeighbors] : distance to each of the closest points
  */
  auto plon = pointLL.lon();
  auto plat = pointLL.lat();
  knnData.m_numNeighbors = 0;
  auto numNeighbors = knnData.maxNeighbors();
  // Initialize distance and index arrays
  knnData.init_indices();
  knnData.init_dist();

  auto const &lons = gps.reg2d->m_lonsReg2d;
  auto const &lats = gps.reg2d->m_latsReg2d;

  auto &params = gps.params;
  long nx = params.dims[0];
  long ny = params.dims[1];
  size_t nxm = params.isCyclic ? nx + 1 : nx;

  if (plon < lons[0]) plon += PI2;
  if (plon > lons[nxm - 1]) plon -= PI2;

  auto pointFound = gps.reg2d->search(plon, plat, knnData, params.searchRadius, params.isCyclic);

  if (pointFound == false && params.extrapolation)
  {
    auto &nbrIndices = knnData.m_indices;
    auto &nbrDistance = knnData.m_dist;
    int searchResult = 0;

    if (numNeighbors < 4)
    {
      size_t nbrIndices4[4];
      double nbrDistance4[4];
      for (size_t i = 0; i < numNeighbors; ++i) nbrIndices4[i] = SIZE_MAX;
      searchResult = grid_search_square_reg2d_NN(nx, ny, nbrIndices4, nbrDistance4, plat, plon, lats, lons);
      if (searchResult < 0)
      {
        for (size_t i = 0; i < numNeighbors; ++i) nbrIndices[i] = nbrIndices4[i];
        for (size_t i = 0; i < numNeighbors; ++i) nbrDistance[i] = nbrDistance4[i];
      }
    }
    else { searchResult = grid_search_square_reg2d_NN(nx, ny, nbrIndices.data(), nbrDistance.data(), plat, plon, lats, lons); }

    if (searchResult < 0) knnData.m_numNeighbors = numNeighbors;
  }
}

void
grid_search_point_unstruct(GridPointsearch &gps, PointLonLat const &pointLL, KnnData &knnData)
{
  /*
    Input variables:

      pointLL : longitude/latitude of the search point

    Output variables:

      knnData.m_indices[numNeighbors] : index of each of the closest points
      knnData.m_dist[numNeighbors] : distance to each of the closest points
  */
  knnData.m_numNeighbors = 0;
  auto numNeighbors = knnData.maxNeighbors();

  // check some more points if distance is the same use the smaller index
  auto numDist = (numNeighbors > 8) ? numNeighbors + 8 : numNeighbors * 2;
  if (numDist > gps.numPoints) numDist = gps.numPoints;

  if (knnData.m_tmpIndices.empty()) knnData.m_tmpIndices.resize(numDist);
  if (knnData.m_tmpDist.empty()) knnData.m_tmpDist.resize(numDist);
  auto &indices = knnData.m_tmpIndices;
  auto &dist = knnData.m_tmpDist;

  auto numIndices = (numNeighbors == 1) ? grid_pointsearch_nearest(gps, pointLL, indices.data(), dist.data())
                                        : grid_pointsearch_qnearest(gps, pointLL, numDist, indices.data(), dist.data());

  numDist = numIndices;
  if (numDist < numNeighbors) numNeighbors = numDist;

  // Initialize distance and index arrays
  knnData.init_indices();
  knnData.init_dist();
  for (size_t i = 0; i < numDist; ++i) knnData.store_distance(indices[i], dist[i], numNeighbors);

  knnData.check_distance();

  if (knnData.m_needCoords)
  {
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), knnData.m_tgtCoord);
    numNeighbors = knnData.m_numNeighbors;
    for (size_t i = 0; i < numNeighbors; ++i)
    {
      gcLLtoXYZ(gps.plons[knnData.m_indices[i]], gps.plats[knnData.m_indices[i]], knnData.m_srcCoords[i]);
    }
  }
}

void
grid_search_point_smooth(GridPointsearch &gps, PointLonLat const &pointLL, KnnData &knnData)
{
  /*
    Input variables:

      pointLL : longitude/latitude of the search point

    Output variables:

      knnData.m_indices[numNeighbors] : index of each of the closest points
      knnData.m_dist[numNeighbors] : distance to each of the closest points
  */
  knnData.m_numNeighbors = 0;
  auto numNeighbors = knnData.maxNeighbors();
  auto checkDistance = (numNeighbors <= 32);

  // check some more points if distance is the same use the smaller index
  auto numDist = checkDistance ? ((numNeighbors > 8) ? numNeighbors + 8 : numNeighbors * 2) : numNeighbors;
  if (numDist > gps.numPoints) numDist = gps.numPoints;

  if (knnData.m_tmpIndices.empty()) knnData.m_tmpIndices.resize(numDist);
  if (knnData.m_tmpDist.empty()) knnData.m_tmpDist.resize(numDist);
  auto &indices = knnData.m_tmpIndices;
  auto &dist = knnData.m_tmpDist;

  auto numIndices = (numNeighbors == 1) ? grid_pointsearch_nearest(gps, pointLL, indices.data(), dist.data())
                                        : grid_pointsearch_qnearest(gps, pointLL, numDist, indices.data(), dist.data());

  numDist = numIndices;

  if (checkDistance)
  {
    if (numDist < numNeighbors) numNeighbors = numDist;

    // Initialize distance and index arrays
    knnData.init_indices(numNeighbors);
    knnData.init_dist(numNeighbors);
    for (size_t i = 0; i < numDist; ++i) knnData.store_distance(indices[i], dist[i], numNeighbors);
  }
  else
  {
    knnData.m_numNeighbors = numDist;
    for (size_t i = 0; i < numDist; ++i) knnData.m_indices[i] = indices[i];
    for (size_t i = 0; i < numDist; ++i) knnData.m_dist[i] = dist[i];
  }

  knnData.check_distance();

  if (knnData.m_needCoords)
  {
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), knnData.m_tgtCoord);
    numNeighbors = knnData.m_numNeighbors;
    for (size_t i = 0; i < numNeighbors; ++i)
    {
      gcLLtoXYZ(gps.plons[knnData.m_indices[i]], gps.plats[knnData.m_indices[i]], knnData.m_srcCoords[i]);
    }
  }
}

void
remap_search_points(RemapSearch &rsearch, PointLonLat const &pointLL, KnnData &knnData)
{
  // clang-format off
  if      (rsearch.srcGrid->type == RemapGridType::Reg2D)   grid_search_point_reg2d(rsearch.gps, pointLL, knnData);
  else if (rsearch.srcGrid->type == RemapGridType::HealPix) grid_search_point_healpix(rsearch.gps, pointLL, knnData);
  else                                                      grid_search_point_unstruct(rsearch.gps, pointLL, knnData);
  // clang-format on
}

static int
quad_cross_products(double plon, double plat, double (&lons)[4], const double (&lats)[4])
{
  int n;

  // clang-format off
  // For consistency, we must make sure all lons are in same 2pi interval
  auto vec1_lon = lons[0] - plon;
  if      (vec1_lon >  PI) lons[0] -= PI2;
  else if (vec1_lon < -PI) lons[0] += PI2;

  for (n = 1; n < 4; ++n)
  {
    vec1_lon = lons[n] - lons[0];
    if      (vec1_lon >  PI) lons[n] -= PI2;
    else if (vec1_lon < -PI) lons[n] += PI2;
  }

  constexpr double crossEps = 1.e-20;
  int scross[4], scrossLast = 0;
  // corner_loop
  for (n = 0; n < 4; ++n)
  {
    int next_n = (n + 1) % 4;
    // Here we take the cross product of the vector making up each box side
    // with the vector formed by the vertex and search point.
    // If all the cross products are positive, the point is contained in the box.
    auto vec1_lat = lats[next_n] - lats[n];
    vec1_lon = lons[next_n] - lons[n];
    auto vec2_lat = plat - lats[n];
    auto vec2_lon = plon - lons[n];

    // Check for 0,2pi crossings
    if      (vec1_lon >  3.0 * PIH) vec1_lon -= PI2;
    else if (vec1_lon < -3.0 * PIH) vec1_lon += PI2;

    if      (vec2_lon >  3.0 * PIH) vec2_lon -= PI2;
    else if (vec2_lon < -3.0 * PIH) vec2_lon += PI2;

    auto crossProduct = vec1_lon * vec2_lat - vec2_lon * vec1_lat;

    // If cross product is less than ZERO, this cell doesn't work
    // 2008-10-16 Uwe Schulzweida: bug fix for cross_product eq zero
    // 2022-06-11 Uwe Schulzweida: replaced zero by crossEps
    scross[n] = (crossProduct < -crossEps) ? -1 : (crossProduct > crossEps) ? 1 : 0;
    if (n == 0) scrossLast = scross[n];
    if ((scross[n] < 0 && scrossLast > 0) || (scross[n] > 0 && scrossLast < 0)) break;
    scrossLast = scross[n];
  }

  if (n >= 4)
  {
    n = 0;
    if      (scross[0] >= 0 && scross[1] >= 0 && scross[2] >= 0 && scross[3] >= 0) n = 4;
    else if (scross[0] <= 0 && scross[1] <= 0 && scross[2] <= 0 && scross[3] <= 0) n = 4;
  }
  // clang-format on

  return n;
}

bool
point_in_quad(bool isCyclic, size_t nx, size_t ny, size_t i, size_t j, SquareCorners &squareCorners, double plon, double plat,
              const double *centerLons, const double *centerLats)
{
  bool searchResult = false;
  size_t ip1 = (i < (nx - 1)) ? i + 1 : isCyclic ? 0 : i;
  size_t jp1 = (j < (ny - 1)) ? j + 1 : j;

  if (i == ip1 || j == jp1) return searchResult;

  size_t idx[4];
  idx[0] = j * nx + i;
  idx[1] = j * nx + ip1;    // east
  idx[2] = jp1 * nx + ip1;  // north-east
  idx[3] = jp1 * nx + i;    // north

  for (int k = 0; k < 4; ++k) squareCorners.lons[k] = centerLons[idx[k]];
  for (int k = 0; k < 4; ++k) squareCorners.lats[k] = centerLats[idx[k]];

  int n = quad_cross_products(plon, plat, squareCorners.lons, squareCorners.lats);

  // If cross products all same sign, we found the location
  if (n >= 4)
  {
    for (int k = 0; k < 4; ++k) squareCorners.indices[k] = idx[k];
    searchResult = true;
  }

  return searchResult;
}

static int
grid_search_square_curv2d(GridPointsearch &gps, RemapGrid *rgrid, SquareCorners &squareCorners, PointLonLat const &pointLL)
{
  /*
    Input variables:

      pointLL : longitude/latitude of the search point

    Output variables:

      srcIndices[4] :  index of each corner point enclosing P
      srcLats[4]    :  latitudes  of the four corner points
      srcLons[4]    :  longitudes of the four corner points
  */
  int searchResult = 0;

  for (int i = 0; i < 4; ++i) squareCorners.indices[i] = 0;

  double dist = 0.0;
  size_t index = 0;
  auto numIndices = grid_pointsearch_nearest(gps, pointLL, &index, &dist);
  if (numIndices > 0)
  {
    auto nx = rgrid->dims[0];
    auto ny = rgrid->dims[1];

    for (int k = 0; k < 4; ++k)
    {
      // Determine neighbor index
      auto j = index / nx;
      auto i = index - j * nx;
      if (k == 0 || k == 2) i = (i > 0) ? i - 1 : rgrid->isCyclic ? nx - 1 : 0;
      if (k == 0 || k == 1) j = (j > 0) ? j - 1 : 0;
      auto plon = pointLL.lon();
      auto plat = pointLL.lat();
      if (point_in_quad(rgrid->isCyclic, nx, ny, i, j, squareCorners, plon, plat, rgrid->centerLons.data(),
                        rgrid->centerLats.data()))
      {
        searchResult = 1;
        return searchResult;
      }
    }
  }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside the grid.
    Fall back to a distance-weighted average of the four closest points. Go ahead and compute weights here,
    but store in srcLats and return -index to prevent the parent routine from computing bilinear weights.
  */
  if (!rgrid->doExtrapolate) return searchResult;

  constexpr size_t numDist = 4;
  numIndices = grid_pointsearch_qnearest(gps, pointLL, numDist, squareCorners.indices, squareCorners.lats);
  if (numIndices == 4)
  {
    for (int i = 0; i < 4; ++i) squareCorners.lats[i] = 1.0 / (squareCorners.lats[i] + TINY);
    double distance = 0.0;
    for (int i = 0; i < 4; ++i) distance += squareCorners.lats[i];
    for (int i = 0; i < 4; ++i) squareCorners.lats[i] /= distance;
    searchResult = -1;
  }

  return searchResult;
}

int
remap_search_square(RemapSearch &rsearch, PointLonLat const &pointLL, SquareCorners &squareCorners)
{
  if (rsearch.srcGrid->type == RemapGridType::Reg2D)
    return grid_search_square_reg2d(rsearch.srcGrid, squareCorners, pointLL.lat(), pointLL.lon());
  else
    return grid_search_square_curv2d(rsearch.gps, rsearch.srcGrid, squareCorners, pointLL);
}
