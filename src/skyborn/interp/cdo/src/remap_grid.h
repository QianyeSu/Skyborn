/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef REMAP_GRID_H
#define REMAP_GRID_H

#include <numbers>

#include "varray.h"
#include "point.h"
#include "mpim_grid/grid_healpix.h"

enum struct RemapGridType
{
  Undefined,
  HealPix,
  Reg2D,
  Unstruct
};

// clang-format off
struct  // RemapGrid
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
RemapGrid
// clang-format on
{
  std::string name{};
  RemapGridType type{ RemapGridType::Undefined };
  int gridID{ -1 };
  int tmpgridID{ -1 };
  int rank{ 0 };           // rank of the grid
  size_t size{ 0 };        // total points on the grid
  size_t numCorners{ 0 };  // number of corners for each grid cell

  int nside{ 0 };
  HpOrder order{ HpOrder::Undef };
  HpParams hpParams{};

  bool needCellCorners{ false };
  bool useCellCorners{ false };  // use corners for bounding boxes

  bool doExtrapolate{ false };
  bool isCyclic{ false };

  size_t dims[2]{};  // size of grid dimension

  int nvgp{ 0 };       // size of vgpm
  Varray<int> vgpm{};  // flag which cells are valid
  Vmask mask{};        // flag which cells participate

  Varray<double> centerLonsReg2d{};  // reg2d lon/lat coordinates for
  Varray<double> centerLatsReg2d{};  // each grid center in radians
  Varray<double> cornerLonsReg2d{};  // reg2d lon/lat coordinates for
  Varray<double> cornerLatsReg2d{};  // each grid corner in radians

  Varray<double> centerLons{};  // lon/lat coordinates for
  Varray<double> centerLats{};  // each grid center in radians
  Varray<double> cornerLons{};  // lon/lat coordinates for
  Varray<double> cornerLats{};  // each grid corner in radians

  Varray<double> cellArea{};  // total area of each grid cell
  Varray<double> cellFrac{};  // fractional area of grid cells participating in remapping

  PointLonLat
  get_lonlat(size_t index) const
  {
    double lon{}, lat{};

    if (this->type == RemapGridType::Reg2D)
    {
      auto nx = this->dims[0];
      auto iy = index / nx;
      auto ix = index - iy * nx;
      lat = this->centerLatsReg2d[iy];
      lon = this->centerLonsReg2d[ix];
      if (lon < 0) lon += 2.0 * std::numbers::pi;
    }
    else if (this->type == RemapGridType::HealPix)
    {
      hp_index_to_lonlat(this->hpParams, index, &lon, &lat);
      // if (lon < 0) lon += 2.0 * std::numbers::pi;
    }
    else
    {
      lat = this->centerLats[index];
      lon = this->centerLons[index];
    }

    return PointLonLat(lon, lat);
  }
};

#endif /* REMAP_GRID_H */
