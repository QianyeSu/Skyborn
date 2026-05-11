/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "remap.h"

static void
nbr_indices_and_distance(int &searchResult, size_t index, double &distMin, double distance, size_t *nbrIndices, double *nbrDistance)
{
  if (distance < distMin)
  {
    for (int k = 0; k < 4; ++k)
    {
      if (distance < nbrDistance[k])
      {
        for (int i = 3; i > k; --i) nbrIndices[i] = nbrIndices[i - 1];
        for (int i = 3; i > k; --i) nbrDistance[i] = nbrDistance[i - 1];
        searchResult = -1;
        nbrIndices[k] = index;
        nbrDistance[k] = distance;
        distMin = nbrDistance[3];
        break;
      }
    }
  }
}

int
grid_search_square_reg2d_NN(size_t nx, size_t ny, size_t *nbrIndices, double *nbrDistance, double plat, double plon,
                            Varray<double> const &lats, Varray<double> const &lons)
{
  // nbrIndices and nbrDistance must have memory for four values!
  int searchResult = 0;

  auto tgtCosLat = std::cos(plat);
  auto tgtSinLat = std::sin(plat);
  auto tgtCosLon = std::cos(plon);
  auto tgtSinLon = std::sin(plon);

  double distMin = DBL_MAX;
  for (int n = 0; n < 4; ++n) nbrDistance[n] = DBL_MAX;

  size_t jjf = 0;
  size_t jjl = ny - 1;
  if (plon >= lons[0] && plon <= lons[nx - 1])
  {
    if (lats[0] < lats[ny - 1])
    {
      if (plat <= lats[0])
        jjl = (ny == 1) ? 0 : 1;
      else
        jjf = (ny == 1) ? 0 : ny - 2;
    }
    else
    {
      if (plat >= lats[0])
        jjl = (ny == 1) ? 0 : 1;
      else
        jjf = (ny == 1) ? 0 : ny - 2;
    }
  }

  std::vector<double> sincosLon(nx);
  for (size_t ii = 0; ii < nx; ++ii) sincosLon[ii] = tgtCosLon * std::cos(lons[ii]) + tgtSinLon * std::sin(lons[ii]);

  for (size_t jj = jjf; jj <= jjl; ++jj)
  {
    auto cosLat = tgtCosLat * std::cos(lats[jj]);
    auto sinLat = tgtSinLat * std::sin(lats[jj]);

    auto jjSkip = (jj > 1 && jj < (ny - 2));

    if (jjSkip)
    {
      size_t ii = 0;
      auto distance = std::acos(cosLat * sincosLon[ii] + sinLat);
      nbr_indices_and_distance(searchResult, jj * nx + ii, distMin, distance, nbrIndices, nbrDistance);
      ii = nx - 1;
      distance = std::acos(cosLat * sincosLon[ii] + sinLat);
      nbr_indices_and_distance(searchResult, jj * nx + ii, distMin, distance, nbrIndices, nbrDistance);
    }
    else
    {
      for (size_t ii = 0; ii < nx; ++ii)
      {
        auto distance = std::acos(cosLat * sincosLon[ii] + sinLat);
        nbr_indices_and_distance(searchResult, jj * nx + ii, distMin, distance, nbrIndices, nbrDistance);
      }
    }
  }

  for (int i = 0; i < 4; ++i) nbrDistance[i] = 1.0 / (nbrDistance[i] + TINY);
  double distance = 0.0;
  for (int i = 0; i < 4; ++i) distance += nbrDistance[i];
  for (int i = 0; i < 4; ++i) nbrDistance[i] /= distance;

  return searchResult;
}

int
grid_search_square_reg2d(RemapGrid *srcGrid, SquareCorners &squareCorners, double plat, double plon)
{
  /*
    Input variables:

      plat : latitude  of the search point
      plon : longitude of the search point

    Output variables:

      indices[4]  : indices of each corner point enclosing P
      lats[4]     : latitudes  of the four corner points
      lons[4]     : longitudes of the four corner points
  */
  int searchResult = 0;
  auto const &latsReg2d = srcGrid->centerLatsReg2d;
  auto const &lonsReg2d = srcGrid->centerLonsReg2d;

  for (int n = 0; n < 4; ++n) squareCorners.indices[n] = 0;

  auto nx = srcGrid->dims[0];
  auto ny = srcGrid->dims[1];

  auto nxm = srcGrid->isCyclic ? nx + 1 : nx;

  if (plon < lonsReg2d[0]) plon += PI2;
  if (plon > lonsReg2d[nxm - 1]) plon -= PI2;

  auto &srcIndices = squareCorners.indices;
  auto &srcLons = squareCorners.lons;
  auto &srcLats = squareCorners.lats;
  size_t ii, jj;
  auto lfound = rect_grid_search(ii, jj, plon, plat, nxm, ny, lonsReg2d, latsReg2d);
  if (lfound)
  {
    size_t iix = (srcGrid->isCyclic && ii == (nxm - 1)) ? 0 : ii;
    srcIndices[0] = (jj - 1) * nx + (ii - 1);
    srcIndices[1] = (jj - 1) * nx + iix;
    srcIndices[2] = jj * nx + iix;
    srcIndices[3] = jj * nx + (ii - 1);

    srcLons[0] = lonsReg2d[ii - 1];
    srcLons[1] = lonsReg2d[iix];
    // For consistency, we must make sure all lons are in same 2pi interval
    if (srcLons[0] > PI2) { srcLons[0] -= PI2; }
    if (srcLons[0] < 0) { srcLons[0] += PI2; }
    if (srcLons[1] > PI2) { srcLons[1] -= PI2; }
    if (srcLons[1] < 0) { srcLons[1] += PI2; }
    srcLons[2] = srcLons[1];
    srcLons[3] = srcLons[0];

    srcLats[0] = latsReg2d[jj - 1];
    srcLats[1] = srcLats[0];
    srcLats[2] = latsReg2d[jj];
    srcLats[3] = srcLats[2];

    searchResult = 1;

    return searchResult;
  }

  /*
    If no cell found, point is likely either in a box that straddles either pole
    or is outside the grid. Fall back to a distance-weighted average of the four
    closest points. Go ahead and compute weights here, but store in srcLats and
    return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!srcGrid->doExtrapolate) return searchResult;

  searchResult = grid_search_square_reg2d_NN(nx, ny, srcIndices, srcLats, plat, plon, latsReg2d, lonsReg2d);

  return searchResult;
}
