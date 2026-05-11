#ifndef GRID_HEALPIX_H
#define GRID_HEALPIX_H

#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <utility>

enum class HpOrder
{
  Undef,
  XY,
  Ring,
  Nested
};

class HpParams
{
public:
  HpParams(){};
  HpParams(int nside, HpOrder order) : m_nside(nside), m_order(order){};

  // clang-format off
  constexpr int nside() const { return m_nside; };
  constexpr int level() const { return static_cast<int>(std::log2(m_nside)); }; // refinement level
  constexpr HpOrder order() const { return m_order; };
  // clang-format on

private:
  int m_nside{ 0 };
  HpOrder m_order{ HpOrder::Undef };
};

HpOrder hp_get_order(std::string const &orderName);
int64_t hp_lonlat_to_index(HpParams hpParams, double xval, double yval);
void hp_index_to_lonlat(HpParams hpParams, int64_t index, double *xval, double *yval);
void hp_get_neighbours(HpParams hpParams, int64_t index, int64_t *neighbours);
void hp_bilinear_interpolate_weights(HpParams hpParams, double lon, double lat, size_t (&indices)[4], double (&weights)[4]);
std::pair<double, double> hp_nested_generate_lonlat(int nside, int64_t index);
void hp_generate_coords(HpOrder order, int nside, int64_t nvals, double *xvals, double *yvals, bool withBounds, double *xbounds,
                        double *ybounds);
void hp_generate_coords_indices(HpOrder order, int nside, int64_t nvals, double *xvals, double *yvals, bool withBounds,
                                double *xbounds, double *ybounds, int64_t *indices);
void hp_generate_latitudes(int nside, std::vector<double> &latitudes);
void hp_generate_ring_indices(HpParams hpParams, size_t gridsize, std::vector<int> &ringIndices, std::vector<int> &ringRows);

template <typename T>
void hp_ring_to_nested(int nside, size_t gridsize, T *arrayIn, T *arrayOut);

template <typename T>
void hp_nested_to_ring(int nside, size_t gridsize, T *arrayIn, T *arrayOut);

#endif
