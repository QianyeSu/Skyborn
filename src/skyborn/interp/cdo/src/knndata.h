/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef KNNDATA_H
#define KNNDATA_H

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <memory>

#include "cdo_math.h"
#include "varray.h"

enum struct WeightingMethod
{
  undefined,
  arithmeticAverage,
  distanceWeighted,
  gaussWeighted,
  rbf,
  linear,
};

std::string weightingMethod_to_string(WeightingMethod method);
WeightingMethod string_to_weightingMethod(std::string const &methodStr);

constexpr double CDO_INTERP_NNN_GAUSS_SCALE_DEFAULT = 0.1;
constexpr double CDO_INTERP_NNN_RBF_SCALE_DEFAULT = 14.87973;

struct KnnParams
{
  size_t k{ 1 };
  size_t kMin{ 0 };
  double maxSearchDistance{ 0.0 };
  double gaussScale{ CDO_INTERP_NNN_GAUSS_SCALE_DEFAULT };
  double rbfScale{ CDO_INTERP_NNN_RBF_SCALE_DEFAULT };
  WeightingMethod weighted{ WeightingMethod::distanceWeighted };
  bool extrapolate{ false };
  // linear
  double searchRadius{ 0.0 };
  double weight0{ 1.0 };
  double weightR{ 1.0 };
};

class KnnData
{
private:
  size_t m_kMin{ 0 };
  size_t m_maxPoints{ 0 };
  size_t m_maxNeighbors{ 0 };
  double m_gaussScale{ CDO_INTERP_NNN_GAUSS_SCALE_DEFAULT };
  double m_rbfScale{ CDO_INTERP_NNN_RBF_SCALE_DEFAULT };
  // linear
  double m_searchRadius{ 0.0 };
  double m_weight0{ 1.0 };
  double m_weightR{ 1.0 };

  void apply_mask(Vmask const &gridMask);

  size_t compute_weights_avg();
  size_t compute_weights_dist();
  size_t compute_weights_linear();
  size_t compute_weights_gauss();
  size_t compute_weights_rbf();

public:
  WeightingMethod m_weighted{ WeightingMethod::distanceWeighted };
  size_t m_numNeighbors{};
  std::vector<size_t> m_indices{};  // source indices at nearest neighbors
  std::vector<double> m_dist{};     // angular distance four nearest neighbors
  std::vector<size_t> m_tmpIndices{};
  std::vector<double> m_tmpDist{};
  std::unique_ptr<double[][3]> m_srcCoords{};
  double m_tgtCoord[3];
  bool m_needCoords{ false };

  inline void
  init()
  {
    m_indices.resize(m_maxNeighbors);
    m_dist.resize(m_maxNeighbors);
    // check some more points if distance is the same use the smaller index
    m_maxPoints = (m_maxNeighbors > 8) ? m_maxNeighbors + 8 : m_maxNeighbors * 2;
    m_needCoords = (m_weighted == WeightingMethod::gaussWeighted || m_weighted == WeightingMethod::rbf);
    if (m_needCoords) { m_srcCoords = std::make_unique<double[][3]>(m_maxPoints); }
  }

  explicit KnnData(KnnParams const &knnParams) : m_weighted(knnParams.weighted)
  {
    m_maxNeighbors = knnParams.k;
    m_kMin = knnParams.kMin;
    m_gaussScale = knnParams.gaussScale;
    m_rbfScale = knnParams.rbfScale;
    m_searchRadius = knnParams.searchRadius;
    m_weight0 = knnParams.weight0;
    m_weightR = knnParams.weightR;

    init();
  }
  explicit KnnData(size_t maxNeighbors) : m_maxNeighbors(maxNeighbors) { init(); }
  explicit KnnData(KnnData &&other) : m_weighted(other.m_weighted)
  {
    m_maxNeighbors = other.m_maxNeighbors;
    m_kMin = other.m_kMin;
    m_gaussScale = other.m_gaussScale;
    m_rbfScale = other.m_rbfScale;
    m_searchRadius = other.m_searchRadius;
    m_weight0 = other.m_weight0;
    m_weightR = other.m_weightR;
    m_needCoords = other.m_needCoords;

    m_maxPoints = other.m_maxPoints;
    m_indices = std::move(other.m_indices);
    m_dist = std::move(other.m_dist);
    m_srcCoords = std::move(other.m_srcCoords);
  }
  ~KnnData() {}

  inline size_t
  maxNeighbors() const
  {
    return m_maxNeighbors;
  }

  inline size_t
  numNeighbors() const
  {
    return m_numNeighbors;
  }

  inline void
  init_indices(size_t numNeighbors)
  {
    for (size_t i = 0; i < numNeighbors; ++i) { m_indices[i] = SIZE_MAX; }
  }

  inline void
  init_dist(size_t numNeighbors)
  {
    for (size_t i = 0; i < numNeighbors; ++i) { m_dist[i] = DBL_MAX; }
  }

  inline void
  init_indices()
  {
    init_indices(m_maxNeighbors);
  }

  inline void
  init_dist()
  {
    init_dist(m_maxNeighbors);
  }

  void store_distance(size_t index, double distance, size_t maxNeighbors);

  inline void
  check_distance()
  {
    constexpr double eps = 1.e-14;
    // If distance is zero, set to small number
    for (size_t i = 0; i < m_numNeighbors; ++i)
      if (m_dist[i] <= 0.0) m_dist[i] = eps;
  }

  size_t
  compute_weights()
  {
    if (m_kMin > 0 && m_numNeighbors < m_kMin) return 0;
    if (m_weighted == WeightingMethod::arithmeticAverage) return compute_weights_avg();
    if (m_weighted == WeightingMethod::distanceWeighted) return compute_weights_dist();
    if (m_weighted == WeightingMethod::linear) return compute_weights_linear();
    if (m_weighted == WeightingMethod::gaussWeighted) return compute_weights_gauss();
    if (m_weighted == WeightingMethod::rbf) return compute_weights_rbf();
    return 0;
  }

  size_t
  compute_weights(Vmask const &gridMask)
  {
    // Compute weights if grid mask is false, eliminate those points
    apply_mask(gridMask);
    return compute_weights();
  }

  template <typename T>
  double
  array_weights_sum(Varray<T> const &array) const
  {
    double sum = 0.0;
    for (size_t i = 0; i < m_numNeighbors; ++i) { sum += array[m_indices[i]] * m_dist[i]; }
    return sum;
  }
};

#endif
