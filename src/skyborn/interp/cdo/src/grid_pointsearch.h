/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef GRID_POINTSEARCH_H
#define GRID_POINTSEARCH_H

#include "pointsearch_reg2d.h"
#include "pointsearch_healpix.h"
#include "pointsearch_unstruct.h"
#include "pointsearch_utils.h"
#include "knndata.h"
#include "remap_grid.h"
#include "cdo_options.h"
#include "varray.h"

#define GPS_NOT_FOUND SIZE_MAX

class GridPointsearch
{
public:
  GridPointsearch() { params.searchRadius = cdo_get_search_radius() * DEG2RAD; }
  ~GridPointsearch()
  {
    if (reg2d) delete reg2d;
    if (healpix) delete healpix;
  }

  void
  set_radius(double chordRadius)
  {
    params.searchRadius = chordRadius;
  }

  void
  enable_extrapolation()
  {
    params.extrapolation = true;
  }

  size_t numPoints = 0;

  const double *plons{ nullptr };
  const double *plats{ nullptr };

  UnstructMethod unstructMethod{ UnstructMethod::nanoflann };

  PointsearchParams params;

  // private:
  PointsearchReg2d *reg2d{ nullptr };
  PointsearchHealpix *healpix{ nullptr };
  PointsearchUnstruct unstruct{};
};

void set_pointsearch_method(std::string const &methodStr);

void grid_search_point_unstruct(GridPointsearch &gps, PointLonLat const &pointLL, KnnData &knnData);
void grid_search_point_smooth(GridPointsearch &gps, PointLonLat const &pointLL, KnnData &knnData);

void grid_pointsearch_create_unstruct(GridPointsearch &gps, Varray<double> const &lons, Varray<double> const &lats,
                                      bool useBoundBox = false);
void grid_pointsearch_create(GridPointsearch &gps, RemapGrid const &remapGrid);
size_t grid_pointsearch_nearest(GridPointsearch &gps, PointLonLat const &pointLL, size_t *index, double *dist);
size_t grid_pointsearch_qnearest(GridPointsearch &gps, PointLonLat const &pointLL, size_t nnn, size_t *indices, double *dist);

struct SquareCorners
{
  double lons[4];     //  longitudes of four corners
  double lats[4];     //  latitudes  of four corners
  size_t indices[4];  //  indices for the four points
};

bool point_in_quad(bool isCyclic, size_t nx, size_t ny, size_t i, size_t j, SquareCorners &squareCorners, double plon, double plat,
                   double const *centerLons, double const *centerLats);

#endif
