/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_NANOFLANN_H
#define POINTSEARCH_NANOFLANN_H

#include "pointsearch_unstruct.h"
#include "pointsearch_utils.h"
#include "cdo_omp.h"
#include "cdo_math.h"
#include "varray.h"
#include "grid_convert.h"
#include "nanoflann.hpp"

template <typename T>
struct PointCloud
{
  struct Point
  {
    T x, y, z;
  };
  std::vector<Point> pts;
  T min[3], max[3];

  // Must return the number of data points
  inline size_t
  kdtree_get_point_count() const
  {
    return pts.size();
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate
  // value, the
  //  "if/else's" are actually solved at compile time.
  inline T
  kdtree_get_pt(size_t idx, int dim) const
  {
    // clang-format off
    if      (dim == 0) return pts[idx].x;
    else if (dim == 1) return pts[idx].y;
    else               return pts[idx].z;
    // clang-format on
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in
  //   "bb" so it can be avoided to redo it again. Look at bb.size() to find out
  //   the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool
  kdtree_get_bbox(BBOX &bb) const
  {
    for (int i = 0; i < 3; ++i) bb[i].low = min[i];
    for (int i = 0; i < 3; ++i) bb[i].high = max[i];
    return true;
  }
  // bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

using nfTree_t
    = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3>;

class PointsearchNanoflann : public PointsearchStrategy
{
public:
  PointsearchNanoflann(Varray<double> const &lons, Varray<double> const &lats, const PointsearchParams &params) : m_params(params)
  {
    create(lons, lats);
  }

  size_t
  search_nearest(PointLonLat const &pointLL, size_t *index, double *dist) override
  {
    if (m_nfTree == nullptr) return 0;

    auto sqrDistMax = cdo::sqr(m_params.searchRadius);

    double tgtPoint[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtPoint);

    if (!m_params.extrapolation)
      for (int i = 0; i < 3; ++i)
        if (tgtPoint[i] < m_min[i] || tgtPoint[i] > m_max[i]) return 0;

    constexpr size_t numResults = 1;
    size_t retIndex;
    double sqrDist;
    nanoflann::KNNResultSet<double> resultSet(sqrDistMax, numResults);
    resultSet.init(&retIndex, &sqrDist);
    m_nfTree->findNeighbors(resultSet, tgtPoint, nanoflann::SearchParams(10));

#define GPS_NOT_FOUND SIZE_MAX
    if (retIndex != GPS_NOT_FOUND)
    {
      *index = retIndex;
      *dist = std::sqrt(sqrDist);
      return 1;
    }

    return 0;
  }

  size_t
  search_qnearest(PointLonLat const &pointLL, size_t nnn, size_t *indices, double *dist) override
  {
    size_t numIndices{ 0 };

    if (m_nfTree == nullptr) return numIndices;

    auto sqrDistMax = cdo::sqr(m_params.searchRadius);

    double tgtPoint[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtPoint);

    if (!m_params.extrapolation)
      for (int i = 0; i < 3; ++i)
        if (tgtPoint[i] < m_min[i] || tgtPoint[i] > m_max[i]) return numIndices;

    numIndices = m_nfTree->knnRangeSearch(&tgtPoint[0], sqrDistMax, nnn, &indices[0], &dist[0]);
    for (size_t i = 0; i < numIndices; ++i) dist[i] = std::sqrt(dist[i]);

    return numIndices;
  }

private:
  float m_min[3]{};
  float m_max[3]{};
  std::unique_ptr<PointCloud<double>> m_pointCloud;
  std::unique_ptr<nfTree_t> m_nfTree;
  const PointsearchParams &m_params;

  void
  create(Varray<double> const &lons, Varray<double> const &lats)
  {
    auto n = lons.size();
    m_pointCloud = std::make_unique<PointCloud<double>>();

    double min[3] = { 1.e9, 1.e9, 1.e9 };
    double max[3] = { -1.e9, -1.e9, -1.e9 };

    // Generating Point Cloud
    m_pointCloud->pts.resize(n);
#ifdef HAVE_OPENMP45
#pragma omp parallel for if (n > cdoMinLoopSize) schedule(static) reduction(min : min[ : 3]) reduction(max : max[ : 3])
#endif
    for (size_t i = 0; i < n; ++i)
    {
      double pointXYZ[3];
      gcLLtoXYZ(lons[i], lats[i], pointXYZ);
      m_pointCloud->pts[i].x = pointXYZ[0];
      m_pointCloud->pts[i].y = pointXYZ[1];
      m_pointCloud->pts[i].z = pointXYZ[2];
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

    for (int i = 0; i < 3; ++i) m_pointCloud->min[i] = min[i];
    for (int i = 0; i < 3; ++i) m_pointCloud->max[i] = max[i];

    // construct a kd-tree index:
    m_nfTree = std::make_unique<nfTree_t>(3 /*dim*/, *m_pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */));
    m_nfTree->buildIndex();
  }
};

#endif
