/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef MPIM_GRID_H
#define MPIM_GRID_H

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

#include <cdi.h>

#include "grid_proj.h"
#include "grid_rot.h"
#include "grid_convert.h"
#include "varray.h"

enum class LonLatUnits
{
  Deg,
  Rad
};

struct CoordRange
{
  double minVal{};
  double maxVal{};
  bool hasNan{};
};

enum class NeedCorners
{
  No = 0,
  Yes = 1,
  IfAvail = 2
};

extern bool gridVerbose;

static inline bool
gridProjIsSupported(int gridID)
{
  auto gridtype = gridInqType(gridID);
  auto pt = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;
  return (pt == CDI_PROJ_RLL || pt == CDI_PROJ_LCC || pt == CDI_PROJ_LAEA || pt == CDI_PROJ_STERE || pt == CDI_PROJ_SINU
          || pt == CDI_PROJ_HEALPIX);
}

static inline bool
gridHasCoordinates(int gridID)
{
  return (gridInqXvals(gridID, nullptr) && gridInqYvals(gridID, nullptr));
}

static inline bool
gridHasBounds(int gridID)
{
  return (gridInqXbounds(gridID, nullptr) && gridInqYbounds(gridID, nullptr));
}

void gridEnableVerbose(bool enable);

int nfc_to_nlat(int nfc, int ntr);
int nlat_to_ntr(int nlat);
int nlat_to_ntr_linear(int nlat);
int nlat_to_ntr_cubic(int nlat);
int ntr_to_nlat(int ntr);
int ntr_to_nlat_linear(int ntr);
int ntr_to_nlat_cubic(int ntr);
int nlat_to_nlon(int nlat);
int nlat_to_nlon_cubic(int nlat);

void grid_copy_names(int gridID1, int gridID2);
void grid_copy_mapping(int gridID1, int gridID2);

bool grid_is_distance_generic(int gridID);

void grid_to_radian(Varray<double> &values);
void cdo_grid_to_radian(int gridID, int varID, Varray<double> &values, std::string const &description);
void cdo_grid_to_degree(int gridID, int varID, Varray<double> &values, std::string const &description);
void cdo_grid_to_degree(int gridID, int varID, size_t nvals, double *values, std::string const &description);

void grid_gen_corners(size_t n, Varray<double> const &vals, Varray<double> &corners);
void grid_gen_bounds(size_t n, Varray<double> const &vals, Varray<double> &bounds);
void grid_check_lat_borders(int n, Varray<double> &ybounds);

void grid_gen_xbounds2D(size_t nx, size_t ny, Varray<double> const &xbounds, Varray<double> &xbounds2D);
void grid_gen_ybounds2D(size_t nx, size_t ny, Varray<double> const &ybounds, Varray<double> &ybounds2D);

int gridcell_weights(int gridID, Varray<double> &weights);
int gridGenArea(int gridID, Varray<double> &area);
int gridGenAreaReg2Dweights(int gridID, Varray<double> &area);

int gridToZonal(int gridID);
int gridToMeridional(int gridID);
int gridToUnstructured(int gridID, NeedCorners needCorners = NeedCorners::No);
int gridToUnstructuredSelecton(int gridID1, std::vector<size_t> const &selectionIndexList, int nocoords, int nobounds);
int gridToCurvilinear(int gridID, NeedCorners needCorners = NeedCorners::No);
int gridCurvilinearToRegular(int gridID);
int gridProjectionToRegular(int gridID);
int gridToRegular(int gridID);
void field2regular(int gridID1, int gridID2, double missval, Varray<double> &array, size_t numMissVals, int lnearest);

// GME grid
void gme_factorni(int kni, int *kni2, int *kni3);
void gme_grid(int withBounds, size_t gridsize, double *rlon, double *rlat, double *blon, double *blat, int *imask, int ni, int nd,
              int ni2, int ni3);

void cdo_print_griddes(int gridID, int opt, bool genBounds = false);

bool grid_has_proj_params(int gridID);
std::string grid_get_proj_params(int gridID);

bool is_point_grid(int gridID);
int generate_full_point_grid(int gridID);
int generate_full_cell_grid(int gridID);
int generate_full_grid(int gridID);

static inline bool
is_unstruct_grid(int gridID)
{
  auto gridType = gridInqType(gridID);
  return (gridType == GRID_UNSTRUCTURED || gridType == GRID_GME || gridType == GRID_HEALPIX
          || (gridType == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_HEALPIX));
}

static inline bool
is_reg2d_grid(int gridID)
{
  auto gridType = gridInqType(gridID);
  return (gridType == GRID_LONLAT || gridType == GRID_GAUSSIAN);
}

static inline bool
is_global_healpix_grid(int gridID)
{
  bool isGlobalHealpix{ false };
  auto gridType = gridInqType(gridID);
  if (gridType == GRID_HEALPIX)
  {
    auto numIndices = gridInqIndices(gridID, nullptr);
    if (numIndices == 0 /*|| numIndices == gridInqSize(gridID)*/) isGlobalHealpix = true;
  }
  return (isGlobalHealpix || (gridType == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_HEALPIX));
}

static inline bool
is_healpix_grid(int gridID)
{
  auto gridType = gridInqType(gridID);
  return (gridType == GRID_HEALPIX || (gridType == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_HEALPIX));
}

static inline double
orthodrome(double x1, double y1, double x2, double y2)
{
  return std::acos(std::sin(y1) * std::sin(y2) + std::cos(y1) * std::cos(y2) * std::cos(x2 - x1));
}

LonLatUnits string_to_LonLatUnits(std::string const &units, std::string const &description);
LonLatUnits cdo_grid_get_units(int gridID, int varID, std::string const &description);

void check_longitude_range(Varray<double> const &xvals, std::string const &txt, LonLatUnits units);
void check_latitude_range(Varray<double> const &yvals, std::string const &txt, LonLatUnits units);

#endif /* MPIM_GRID_H */
