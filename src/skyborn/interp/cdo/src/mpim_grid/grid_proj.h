/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef GRID_PROJ_H
#define GRID_PROJ_H

#include <cstddef>
#include <string>
#include "varray.h"

std::vector<size_t> healpix_compute_cell_parameter(int refinementLevel, double lon1, double lon2, double lat1, double lat2);

int cdo_healpix_to_lonlat(int gridID, size_t nvals, Varray<double> &xvals, Varray<double> &yvals, bool withBounds,
                          Varray<double> &xbounds, Varray<double> &ybounds);

int cdo_lcc_to_lonlat(int gridID, size_t nvals, Varray<double> &xvals, Varray<double> &yvals);
int cdo_stere_to_lonlat(int gridID, size_t nvals, Varray<double> &xvals, Varray<double> &yvals);

void cdo_sinu_to_lonlat(size_t nvals, Varray<double> &xvals, Varray<double> &yvals);
void cdo_laea_to_lonlat(int gridID, size_t nvals, Varray<double> &xvals, Varray<double> &yvals);

void cdo_proj_to_lonlat(std::string const &projParams, size_t nvals, Varray<double> &xvals, Varray<double> &yvals);

int proj_lonlat_to_lcc(struct CDI_GridProjParams gpp, size_t nvals, double *xvals, double *yvals);
int proj_lcc_to_lonlat(struct CDI_GridProjParams gpp, double x_0, double y_0, size_t nvals, double *xvals, double *yvals);

int proj_lonlat_to_stere(struct CDI_GridProjParams gpp, size_t nvals, double *xvals, double *yvals);
int proj_stere_to_lonlat(struct CDI_GridProjParams gpp, double x_0, double y_0, size_t nvals, double *xvals, double *yvals);

void grid_def_params_laea(int gridID, double a, double lon0, double lat0);
void grid_def_params_sinu(int gridID);

#endif /* GRID_PROJ_H */
