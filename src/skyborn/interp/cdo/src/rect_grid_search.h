/*
  Narrow helper extracted from the upstream CDO interpol search utilities.
  Keep this local header limited to the regular-grid interval lookup helpers
  that are directly needed by the vendored remap subset.
*/

#ifndef SKYBORN_CDO_RECT_GRID_SEARCH_H
#define SKYBORN_CDO_RECT_GRID_SEARCH_H

#include <cstddef>
#include "varray.h"

bool rect_grid_search(
    size_t &ii,
    size_t &jj,
    double x,
    double y,
    size_t nxm,
    size_t nym,
    Varray<double> const &xm,
    Varray<double> const &ym
);

bool rect_grid_search2(
    long &imin,
    long &imax,
    double xmin,
    double xmax,
    long nxm,
    Varray<double> const &xm
);

#endif
