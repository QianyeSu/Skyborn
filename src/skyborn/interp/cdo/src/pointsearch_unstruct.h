/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_H
#define POINTSEARCH_H

#include <cstdio>
#include <memory>
#include "point.h"

enum struct UnstructMethod
{
  undefined,
  kdtree,
  nanoflann,
  spherepart,
  full,
};

class PointsearchStrategy
{
public:
  virtual ~PointsearchStrategy() = default;
  virtual size_t search_nearest(PointLonLat const &pointLL, size_t *index, double *dist) = 0;
  virtual size_t search_qnearest(PointLonLat const &pointLL, size_t nnn, size_t *index, double *dist) = 0;
};

class PointsearchUnstruct
{
public:
  // explicit PointsearchUnstruct(std::unique_ptr<PointsearchStrategy> &&strategy = {}) : m_strategy(std::move(strategy)) {}
  // PointsearchUnstruct(std::unique_ptr<PointsearchStrategy> &&strategy = {}) : m_strategy(std::move(strategy)) {}
  // PointsearchUnstruct() {}
  ~PointsearchUnstruct()
  {
    if (m_strategy) delete m_strategy;
  }
  /*
  void
  set_strategy(std::unique_ptr<PointsearchStrategy> &&strategy)
  {
    m_strategy = std::move(strategy);
  }
  */
  void
  set_strategy(PointsearchStrategy *strategy)
  {
    if (m_strategy) delete m_strategy;
    m_strategy = strategy;
  }

  size_t
  search_nearest(PointLonLat const &pointLL, size_t *index, double *dist)
  {
    if (m_strategy) { return m_strategy->search_nearest(pointLL, index, dist); }
    std::fprintf(stderr, "PointsearchUnstruct::search_nearest: PointsearchStrategy not initialized!\n");
    return 0;
  }

  size_t
  search_qnearest(PointLonLat const &pointLL, size_t nnn, size_t *index, double *dist)
  {
    if (m_strategy) { return m_strategy->search_qnearest(pointLL, nnn, index, dist); }
    std::fprintf(stderr, "PointsearchUnstruct::search_qnearest: PointsearchStrategy not initialized!\n");
    return 0;
  }

private:
  // std::unique_ptr<PointsearchStrategy> m_strategy{};
  PointsearchStrategy *m_strategy{ nullptr };
};

#endif
