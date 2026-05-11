/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_SPHEREPART_H
#define POINTSEARCH_SPHEREPART_H

#include "pointsearch_unstruct.h"
#include "pointsearch_utils.h"
#include "varray.h"
#include "cdo_omp.h"
#include "grid_convert.h"
extern "C"
{
#include "lib/yac/src/grid_cell.h"
#include "lib/yac/src/sphere_part.h"
}

class PointsearchSpherepart : public PointsearchStrategy
{
public:
  PointsearchSpherepart(Varray<double> const &lons, Varray<double> const &lats, const PointsearchParams &params) : m_params(params)
  {
    create(lons, lats);
  }
  ~PointsearchSpherepart() override
  {
    if (m_yacPointSearch) yac_delete_point_sphere_part_search(m_yacPointSearch);
  }

  size_t
  search_nearest(PointLonLat const &pointLL, size_t *index, double *dist) override
  {
    if (m_yacPointSearch == nullptr) return 0;

    auto searchArcRadius = chord_to_arc_length(m_params.searchRadius);

    double tgtPoint[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtPoint);

    if (!m_params.extrapolation)
      for (int i = 0; i < 3; ++i)
        if (tgtPoint[i] < m_min[i] || tgtPoint[i] > m_max[i]) return 0;

    size_t local_point_ids_array_size = 0;
    size_t num_local_point_ids;
    size_t *local_point_ids = nullptr;
    double cos_angle;

    yac_point_sphere_part_search_NN(m_yacPointSearch, 1, &tgtPoint, &cos_angle, nullptr, nullptr, &local_point_ids,
                                    &local_point_ids_array_size, &num_local_point_ids);

    size_t numIndices = 0;
    if (num_local_point_ids > 0)
    {
      *dist = std::acos(cos_angle);
      if (*dist <= searchArcRadius)
      {
        numIndices = 1;
        *index = local_point_ids[0];
        for (size_t i = 1; i < num_local_point_ids; ++i)
          if (local_point_ids[i] < *index) *index = local_point_ids[i];
      }
    }

    if (local_point_ids) free(local_point_ids);

    return numIndices;
  }

  size_t
  search_qnearest(PointLonLat const &pointLL, size_t nnn, size_t *indices, double *dist) override
  {
    size_t numIndices = 0;

    if (m_yacPointSearch == nullptr) return numIndices;

    auto searchArcRadius = chord_to_arc_length(m_params.searchRadius);

    double tgtPoint[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtPoint);

    if (!m_params.extrapolation)
      for (int i = 0; i < 3; ++i)
        if (tgtPoint[i] < m_min[i] || tgtPoint[i] > m_max[i]) return numIndices;

    size_t local_point_ids_array_size = 0;
    size_t num_local_point_ids;
    size_t *local_point_ids = nullptr;

    size_t cos_angles_array_size = 0;
    double *cos_angles = nullptr;

    yac_point_sphere_part_search_NNN(m_yacPointSearch, 1, &tgtPoint, nnn, &cos_angles, &cos_angles_array_size, nullptr, nullptr,
                                     &local_point_ids, &local_point_ids_array_size, &num_local_point_ids);

    if (num_local_point_ids > 0)
    {
      auto maxIndices = (num_local_point_ids < nnn) ? num_local_point_ids : nnn;
      numIndices = 0;
      for (size_t i = 0; i < maxIndices; ++i)
      {
        auto angle = std::acos(cos_angles[i]);
        if (angle < searchArcRadius)
        {
          indices[numIndices] = local_point_ids[i];
          dist[numIndices] = angle;
          numIndices++;
        }
      }
    }

    if (cos_angles) free(cos_angles);
    if (local_point_ids) free(local_point_ids);

    return numIndices;
  }

private:
  float m_min[3] = { 0 };
  float m_max[3] = { 0 };
  std::unique_ptr<double[][3]> m_pointsXYZ;

  point_sphere_part_search *m_yacPointSearch{ nullptr };
  const PointsearchParams &m_params;

  void
  create(Varray<double> const &lons, Varray<double> const &lats)
  {
    auto n = lons.size();
    m_pointsXYZ = std::make_unique<double[][3]>(n);

    double min[3] = { 1.e9, 1.e9, 1.e9 };
    double max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for if (n > cdoMinLoopSize) schedule(static) reduction(min : min[ : 3]) reduction(max : max[ : 3])
#endif
    for (size_t i = 0; i < n; ++i)
    {
      auto &pointXYZ = m_pointsXYZ[i];
      gcLLtoXYZ(lons[i], lats[i], pointXYZ);
      min_point(min, pointXYZ);
      max_point(max, pointXYZ);
    }

    if (!m_params.useBoundBox) min[0] = min[1] = min[2] = -1;
    if (!m_params.useBoundBox) max[0] = max[1] = max[2] = 1;

    adjust_bbox_min(min);
    adjust_bbox_max(max);
    for (int i = 0; i < 3; ++i) m_min[i] = min[i];
    for (int i = 0; i < 3; ++i) m_max[i] = max[i];

    // if (Options::cdoVerbose) cdo_print("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

    auto global_ids = std::make_unique<size_t[]>(n);
    for (size_t i = 0; i < n; ++i) global_ids[i] = i;
    m_yacPointSearch = yac_point_sphere_part_search_new(n, m_pointsXYZ.get(), global_ids.get());
  }
};

#endif
