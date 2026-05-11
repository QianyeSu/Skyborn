/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_HEALPIX_H
#define POINTSEARCH_HEALPIX_H

#include "cdo_math.h"
#include "point.h"
#include "knndata.h"
#include "grid_convert.h"
#include "grid_healpix.h"

class PointsearchHealpix
{
public:
  explicit PointsearchHealpix(const HpParams &hpParams) : m_hpParams(hpParams) {}
  ~PointsearchHealpix() {}

  /* not used because result differ from search() and unstruct::search()
  void
  search_1nn(PointLonLat const &pointLL, knnDataType &knnData, double searchRadius)
  {
    size_t index = (size_t) hp_lonlat_to_index(m_hpParams, pointLL.get_lon(), pointLL.get_lat());

    size_t numIndices = 1;

    double lon, lat;
    hp_index_to_lonlat(m_hpParams, index, &lon, &lat);

    double distance;
    store_distance_healpix(searchRadius, pointLL, knnData, numIndices, &index, &distance, &lon, &lat);
   }
   */

  void
  search(PointLonLat const &pointLL, KnnData &knnData, double searchRadius)
  {
    auto index = hp_lonlat_to_index(m_hpParams, pointLL.lon(), pointLL.lat());

    int64_t neighbours[8];
    hp_get_neighbours(m_hpParams, index, neighbours);

    size_t indices[9];
    indices[0] = index;
    size_t numIndices = 1;
    for (int i = 0; i < 8; ++i)
      if (neighbours[i] >= 0) indices[numIndices++] = neighbours[i];

    double lons[9], lats[9];
    for (size_t i = 0; i < numIndices; ++i) hp_index_to_lonlat(m_hpParams, indices[i], &lons[i], &lats[i]);

    store_distance_healpix(searchRadius, pointLL, knnData, numIndices, indices, lons, lats);

    if (knnData.m_needCoords)
    {
      gcLLtoXYZ(pointLL.lon(), pointLL.lat(), knnData.m_tgtCoord);
      auto numNeighbors = knnData.m_numNeighbors;
      for (size_t i = 0; i < numNeighbors; ++i)
      {
        double lon, lat;
        hp_index_to_lonlat(m_hpParams, knnData.m_indices[i], &lon, &lat);
        gcLLtoXYZ(lon, lat, knnData.m_srcCoords[i]);
      }
    }
  }

private:
  HpParams m_hpParams;

  void
  store_distance_healpix(double searchRadius, PointLonLat const &pointLL, KnnData &knnData, size_t numIndices, size_t *indices,
                         double const *lons, double const *lats)
  {
    double tgtCoord[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtCoord);

    auto sqrSearchRadius = cdo::sqr(searchRadius);

    double distances[9];
    double srcCoord[3];
    size_t numWeights = 0;
    for (size_t i = 0; i < numIndices; ++i)
    {
      gcLLtoXYZ(lons[i], lats[i], srcCoord);
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
  }
};

#endif
