/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_KDTREE_H
#define POINTSEARCH_KDTREE_H

#include "pointsearch_unstruct.h"
#include "pointsearch_utils.h"
#include "cdo_omp.h"
#include "cdo_options.h"
#include "cdo_math.h"
#include "varray.h"
#include "grid_convert.h"
#include "kdtreelib/kdtree.h"

class PointsearchKDtree : public PointsearchStrategy
{
public:
  PointsearchKDtree(Varray<double> const &lons, Varray<double> const &lats, const PointsearchParams &params) : m_params{ params }
  {
    create(lons, lats);
  }
  ~PointsearchKDtree() override
  {
    if (m_kdtree) kd_destroyTree(m_kdtree);
  }

  size_t
  search_nearest(PointLonLat const &pointLL, size_t *index, double *dist) override
  {
    if (m_kdtree == nullptr) return 0;

    auto sqrDistMax = cdo::sqr(m_params.searchRadius);
    auto sqrDist = sqrDistMax;

    double tgtPoint[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtPoint);

    if (!m_params.extrapolation)
      for (int i = 0; i < 3; ++i)
        if (tgtPoint[i] < m_min[i] || tgtPoint[i] > m_max[i]) return 0;

    const auto node = kd_nearest(m_kdtree->node, tgtPoint, &sqrDist, 3);
    if (node && sqrDist < sqrDistMax)
    {
      *index = node->index;
      *dist = std::sqrt(sqrDist);
      return 1;
    }

    return 0;
  }

  size_t
  search_qnearest(PointLonLat const &pointLL, size_t nnn, size_t *indices, double *dist) override
  {
    size_t numIndices = 0;

    if (m_kdtree == nullptr) return numIndices;

    auto sqrDistMax = cdo::sqr(m_params.searchRadius);

    kdata_t tgtPoint[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtPoint);

    if (!m_params.extrapolation)
      for (int i = 0; i < 3; ++i)
        if (tgtPoint[i] < m_min[i] || tgtPoint[i] > m_max[i]) return numIndices;

    kdata_t sqrDist = sqrDistMax;
    auto result = kd_qnearest(m_kdtree->node, tgtPoint, &sqrDist, nnn, 3);
    if (result)
    {
      resItem *p;
      while (pqremove_min(result, &p))
      {
        if (p->dist_sq < sqrDistMax)
        {
          indices[numIndices] = p->node->index;
          dist[numIndices] = std::sqrt(p->dist_sq);
          numIndices++;
        }

        free(p);  // Free the result node taken from the heap
      }
      free(result->d);  // free the heap
      free(result);     // and free the heap information structure
    }

    return numIndices;
  }

private:
  float m_min[3]{};
  float m_max[3]{};
  kdTree_t *m_kdtree{ nullptr };
  const PointsearchParams &m_params;

  void
  create(Varray<double> const &lons, Varray<double> const &lats)
  {
    auto n = lons.size();
    std::vector<kd_point> pointlist(n);
    // see  example_cartesian.c

    kdata_t min[3] = { 1.e9, 1.e9, 1.e9 };
    kdata_t max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for if (n > cdoMinLoopSize) schedule(static) reduction(min : min[ : 3]) reduction(max : max[ : 3])
#endif
    for (size_t i = 0; i < n; ++i)
    {
      auto &point = pointlist[i].point;
      gcLLtoXYZ(lons[i], lats[i], point);
      min_point(min, point);
      max_point(max, point);
      pointlist[i].index = i;
    }

    if (!m_params.useBoundBox) min[0] = min[1] = min[2] = -1;
    if (!m_params.useBoundBox) max[0] = max[1] = max[2] = 1;

    adjust_bbox_min(min);
    adjust_bbox_max(max);
    for (int i = 0; i < 3; ++i) m_min[i] = min[i];
    for (int i = 0; i < 3; ++i) m_max[i] = max[i];

    // if (Options::cdoVerbose) cdo_print("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

    m_kdtree = kd_buildTree(pointlist.data(), n, min, max, 3, Threading::ompNumMaxThreads);
    // if (m_kdtree == nullptr) cdo_abort("kd_buildTree failed!");
  }
};

#endif
