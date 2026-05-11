/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_math.h"
#include "knndata.h"
#include "pointsearch_reg2d.h"
#include <cmath>

static size_t
fill_src_indices(bool isCyclic, long nx, long ny, long ii, long jj, long k, size_t *psrcIndices)
{
  k /= 2;

  auto j0 = jj - k;
  auto jn = jj + k;
  auto i0 = ii - k;
  auto in = ii + k;
  if (j0 < 0) j0 = 0;
  if (jn >= ny) jn = ny - 1;
  if ((in - i0) > nx)
  {
    i0 = 0;
    in = nx - 1;
  }

  size_t numIndices = 0;

  for (long j = j0; j <= jn; ++j)
    for (long i = i0; i <= in; ++i)
    {
      auto ix = i;
      if (isCyclic && ix < 0) ix += nx;
      if (isCyclic && ix >= nx) ix -= nx;
      if (ix >= 0 && ix < nx && j < ny) psrcIndices[numIndices++] = j * nx + ix;
    }

  return numIndices;
}

void
PointsearchReg2d::compute_point(size_t index, double (&xyz)[3])
{
  auto iy = index / m_nx;
  auto ix = index - iy * m_nx;
  xyz[0] = m_cosLats[iy] * m_cosLons[ix];
  xyz[1] = m_cosLats[iy] * m_sinLons[ix];
  xyz[2] = m_sinLats[iy];
}

void
PointsearchReg2d::store_distance_reg2d(double plon, double plat, KnnData &knnData, size_t numIndices, size_t *indices,
                                       double *distances, double searchRadius)
{
  double tgtCoord[3];
  gcLLtoXYZ(plon, plat, tgtCoord);
  auto sqrSearchRadius = cdo::sqr(searchRadius);

  double srcCoord[3];
  size_t numWeights = 0;
  for (size_t i = 0; i < numIndices; ++i)
  {
    compute_point(indices[i], srcCoord);
    // Find distance to this point
    double sqrDistance = (float) cdo::sqr_distance(tgtCoord, srcCoord);
    // Store the index and distance if this is one of the smallest so far
    if (sqrDistance <= sqrSearchRadius)
    {
      indices[numWeights] = indices[i];
      distances[numWeights] = std::sqrt(sqrDistance);
      numWeights++;
    }
  }

  auto maxNeighbors = std::min(numWeights, knnData.maxNeighbors());
  for (size_t i = 0; i < numWeights; ++i) { knnData.store_distance(indices[i], distances[i], maxNeighbors); }

  knnData.check_distance();

  if (knnData.m_needCoords)
  {
    gcLLtoXYZ(plon, plat, knnData.m_tgtCoord);
    auto numNeighbors = knnData.m_numNeighbors;
    for (size_t i = 0; i < numNeighbors; ++i) { compute_point(knnData.m_indices[i], knnData.m_srcCoords[i]); }
  }
}

bool
PointsearchReg2d::search(double plon, double plat, KnnData &knnData, double searchRadius, bool isCyclic)
{
  auto nxm = isCyclic ? m_nx + 1 : m_nx;
  auto numNeighbors = knnData.maxNeighbors();

  size_t ii, jj;
  auto pointFound = rect_grid_search(ii, jj, plon, plat, nxm, m_ny, m_lonsReg2d, m_latsReg2d);
  if (not pointFound) return pointFound;

  if (isCyclic && ii == (nxm - 1)) ii = 0;

  constexpr size_t MAX_SEARCH_CELLS = 25;
  size_t srcIndices[MAX_SEARCH_CELLS];
  double srcDist[MAX_SEARCH_CELLS];
  auto *psrcIndices = srcIndices;
  auto *psrcDist = srcDist;

  size_t k;
  for (k = 3; k < 10000; k += 2)
    if (numNeighbors <= (size_t) (k - 2) * (k - 2)) break;

  std::vector<size_t> tmpIndices;
  std::vector<double> tmpDist;
  if ((k * k) > MAX_SEARCH_CELLS)
  {
    tmpIndices.resize(k * k);
    tmpDist.resize(k * k);
    psrcIndices = tmpIndices.data();
    psrcDist = tmpDist.data();
  }

  auto numIndices = fill_src_indices(isCyclic, m_nx, m_ny, ii, jj, k, psrcIndices);

  store_distance_reg2d(plon, plat, knnData, numIndices, psrcIndices, psrcDist, searchRadius);

  return pointFound;
}
