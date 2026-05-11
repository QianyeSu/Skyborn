/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_REG2D_H
#define POINTSEARCH_REG2D_H

#include <numbers>
#include "pointsearch_utils.h"
#include "knndata.h"
#include "grid_convert.h"
#include "rect_grid_search.h"
#include "varray.h"

class PointsearchReg2d
{
public:
  PointsearchReg2d(Varray<double> const &lons, Varray<double> const &lats, const PointsearchParams &params)
  {
    m_nx = params.dims[0];
    m_ny = params.dims[1];
    create(lons, lats, params.isCyclic);
  }
  ~PointsearchReg2d() {}

  void
  create(Varray<double> const &lons, Varray<double> const &lats, bool isCyclic)
  {
    auto nxm = isCyclic ? m_nx + 1 : m_nx;

    m_lonsReg2d.resize(nxm);
    m_latsReg2d.resize(m_ny);

    varray_copy(nxm, lons, m_lonsReg2d);
    varray_copy(m_ny, lats, m_latsReg2d);

    m_cosLons.resize(m_nx);
    m_sinLons.resize(m_nx);
    m_cosLats.resize(m_ny);
    m_sinLats.resize(m_ny);

    constexpr double pi2 = 2.0 * std::numbers::pi;
    for (size_t n = 0; n < m_nx; ++n)
      {
        auto lon = lons[n];
        if (lon > pi2) lon -= pi2;
        if (lon < 0) lon += pi2;
        m_cosLons[n] = std::cos(lon);
        m_sinLons[n] = std::sin(lon);
      }

    for (size_t n = 0; n < m_ny; ++n)
      {
        auto lat = lats[n];
        m_cosLats[n] = std::cos(lat);
        m_sinLats[n] = std::sin(lat);
      }
  }

  bool search(double plon, double plat, KnnData &knnData, double searchRadius, bool isCyclic);

  void store_distance_reg2d(double plon, double plat, KnnData &knnData, size_t numIndices, size_t *psrcIndices, double *psrcDist,
                            double searchRadius);

  Varray<double> m_lonsReg2d, m_latsReg2d;

private:
  size_t m_nx{ 0 };
  size_t m_ny{ 0 };
  Varray<double> m_cosLats, m_sinLats;  // cosine, sine of grid lats (for distance)
  Varray<double> m_cosLons, m_sinLons;  // cosine, sine of grid lons (for distance)

  void compute_point(size_t index, double (&xyz)[3]);
};

#endif
