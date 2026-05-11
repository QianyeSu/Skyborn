/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_FULL_H
#define POINTSEARCH_FULL_H

#include "pointsearch_unstruct.h"
#include "pointsearch_utils.h"
#include "cdo_omp.h"
#include "cdo_math.h"
#include "cdo_options.h"
#include "varray.h"
#include "grid_convert.h"
#include "kdtreelib/kdtree.h"

class PointsearchFull : public PointsearchStrategy
{
public:
  PointsearchFull(Varray<double> const &lons, Varray<double> const &lats, const PointsearchParams &params) : m_params{ params }
  {
    create(lons, lats);
  }
  //  ~PointsearchFull() {}

  size_t
  search_nearest(PointLonLat const &pointLL, size_t *index, double *dist) override
  {
    if (m_pointsXYZ == nullptr) return 0;

    auto sqrDistMax = cdo::sqr(m_params.searchRadius);

    double tgtPoint[3];
    gcLLtoXYZ(pointLL.lon(), pointLL.lat(), tgtPoint);

    auto closestPoint = m_n;
    double sqrDist = FLT_MAX;
    for (size_t i = 0; i < m_n; ++i)
    {
      double d = (float) cdo::sqr_distance(tgtPoint, m_pointsXYZ[i]);
      if (closestPoint >= m_n || d < sqrDist || (d <= sqrDist && i < closestPoint))
      {
        sqrDist = d;
        closestPoint = i;
      }
    }

    if (closestPoint < m_n && sqrDist < sqrDistMax)
    {
      *index = closestPoint;
      *dist = std::sqrt(sqrDist);
      return 1;
    }

    return 0;
  }

  size_t
  search_qnearest(PointLonLat const &pointLL, size_t nnn, size_t *indices, double *dist) override
  {
    (void) pointLL;
    (void) nnn;
    (void) indices;
    (void) dist;

    static auto warning{ true };
    if (warning)
    {
      warning = false;
      std::fprintf(stderr, "PointsearchFull::search_qnearest() not implemented\n");
    }

    size_t numIndices = 0;
    if (m_pointsXYZ == nullptr) return numIndices;
    return numIndices;
  }

private:
  size_t m_n{ 0 };
  std::unique_ptr<double[][3]> m_pointsXYZ;
  const PointsearchParams &m_params;

  void
  create(Varray<double> const &lons, Varray<double> const &lats)
  {
    auto n = lons.size();
    m_pointsXYZ = std::make_unique<double[][3]>(n);

#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
    for (size_t i = 0; i < n; ++i) { gcLLtoXYZ(lons[i], lats[i], m_pointsXYZ[i]); }

    m_n = n;
  }
};

#endif
