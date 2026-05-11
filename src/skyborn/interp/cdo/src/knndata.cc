#include <cassert>
#include <cstddef>
#include "knndata.h"
#include "cdo_output.h"
#include "interpol.h"

std::string
weightingMethod_to_string(WeightingMethod method)
{
  if (method == WeightingMethod::arithmeticAverage) return "avg";
  if (method == WeightingMethod::distanceWeighted) return "dist";
  if (method == WeightingMethod::linear) return "linear";
  if (method == WeightingMethod::gaussWeighted) return "gauss";
  if (method == WeightingMethod::rbf) return "rbf";

  return "";
}

WeightingMethod
string_to_weightingMethod(std::string const &methodStr)
{
  if (methodStr == "avg") return WeightingMethod::arithmeticAverage;
  if (methodStr == "dist") return WeightingMethod::distanceWeighted;
  if (methodStr == "linear") return WeightingMethod::linear;
  if (methodStr == "gauss") return WeightingMethod::gaussWeighted;
  if (methodStr == "rbf") return WeightingMethod::rbf;

  cdo_abort("method=%s unsupported (available: avg|dist|linear|gauss|rbf)", methodStr);

  return WeightingMethod::undefined;
}

void
KnnData::apply_mask(Vmask const &gridMask)
{
  size_t n = 0;
  for (size_t i = 0; i < m_numNeighbors; ++i)
  {
    if (gridMask[m_indices[i]])
    {
      m_indices[n] = m_indices[i];
      m_dist[n] = m_dist[i];
      n++;
    }
  }

  m_numNeighbors = n;
}

inline bool
distance_is_less(double distance, double distx, size_t index, size_t index2)
{
  constexpr double cmpTolerance = 1.e-12;
  // return (distance < distx || (distance <= distx && index < index2));
  return (distance + cmpTolerance) < distx || (index < index2 && std::fabs(distance - distx) < cmpTolerance);
}

void
KnnData::store_distance(size_t index, double distance, size_t maxNeighbors)
{
  assert(maxNeighbors <= m_maxNeighbors);
  m_numNeighbors = maxNeighbors;

  if (maxNeighbors == 1)
  {
    if (distance_is_less(distance, m_dist[0], index, m_indices[0]))
    {
      m_indices[0] = index;
      m_dist[0] = distance;
    }
  }
  else
  {
    for (size_t i = 0; i < maxNeighbors; ++i)
    {
      if (distance_is_less(distance, m_dist[i], index, m_indices[i]))
      {
        for (size_t n = maxNeighbors - 1; n > i; --n)
        {
          m_indices[n] = m_indices[n - 1];
          m_dist[n] = m_dist[n - 1];
        }
        m_indices[i] = index;
        m_dist[i] = distance;
        break;
      }
    }
  }
}

size_t
KnnData::compute_weights_avg()
{
  if (m_numNeighbors)
  {
    double weight = 1.0 / m_numNeighbors;
    for (size_t i = 0; i < m_numNeighbors; ++i) { m_dist[i] = weight; }
  }

  return m_numNeighbors;
}

size_t
KnnData::compute_weights_dist()
{
  // Compute weights based on inverse distance
  double distTotal = 0.0;  // sum of neighbor distances (for normalizing)
  for (size_t i = 0; i < m_numNeighbors; ++i)
  {
    m_dist[i] = 1.0 / m_dist[i];
    distTotal += m_dist[i];
  }
  // Normalize weights
  for (size_t i = 0; i < m_numNeighbors; ++i) { m_dist[i] = m_dist[i] / distTotal; }

  return m_numNeighbors;
}

size_t
KnnData::compute_weights_linear()
{
  // Compute weights based on linear interpolation
  double distTotal = 0.0;  // sum of neighbor distances (for normalizing)
  for (size_t i = 0; i < m_numNeighbors; ++i)
  {
    m_dist[i] = intlin(m_dist[i], m_weight0, 0, m_weightR, m_searchRadius);
    distTotal += m_dist[i];
  }
  // Normalize weights
  for (size_t i = 0; i < m_numNeighbors; ++i) { m_dist[i] = m_dist[i] / distTotal; }

  return m_numNeighbors;
}

// code from yac routine compute_weights_gauss()
size_t
KnnData::compute_weights_gauss()
{
  double gauss_scale = m_gaussScale;
  auto n = m_numNeighbors;
  std::vector<double> weights(n);
  for (size_t i = 0; i < n; ++i) { weights[i] = m_dist[i] * m_dist[i]; }

  // a) compute sum of source point distances
  double src_distances_sum = 0.0;
  double src_distances_count = 0.5 * (double) (n * n - n);
  for (size_t i = 0; i < n - 1; ++i)
    for (size_t j = i + 1; j < n; ++j) src_distances_sum += std::sqrt(cdo::sqr_distance(m_srcCoords[i], m_srcCoords[j]));

  // b) C = -1 / (c * d_mean^2)
  double src_distances_mean = src_distances_sum / src_distances_count;
  double scale = -1.0 / (gauss_scale * src_distances_mean * src_distances_mean);

  // c) calculate weights
  // w_i = e^(-d_i^2/(c*s^2))
  // w_i = e^(C * d_i^2)
  double weights_sum = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    weights[i] = std::exp(scale * weights[i]);
    weights_sum += weights[i];
  }

  // If the sum of the weights is very low, which can happen in case
  // the target point is very far away from the group source points.
  if (std::fabs(weights_sum) < 1e-9)
  {
    // Due to limited accuracy the exact contribution of each source
    // point cannot be computed. Therefore, the normalisation would
    // generate NaN's. Hence we fall back to inverse distance weighted
    // averge for this target point.
    compute_weights_dist();
    return n;
  }

  // d) scale weights such that SUM(w_i) == 1
  for (size_t i = 0; i < n; ++i) m_dist[i] = weights[i] / weights_sum;

  return n;
}

extern "C"
{
#include "lib/yac/src/compute_weights.h"
}

size_t
KnnData::compute_weights_rbf()
{
  size_t n = m_numNeighbors;
  yac_compute_weights_rbf(m_tgtCoord, m_srcCoords.get(), n, m_dist.data(), m_rbfScale);
  return n;
}
