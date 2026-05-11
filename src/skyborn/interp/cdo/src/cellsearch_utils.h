/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef CELLSEARCH_UTILS_H
#define CELLSEARCH_UTILS_H

#include <stddef.h>

extern "C"
{
#include "lib/yac/src/grid_cell.h"
}

struct GridCell
{
  double *coordinatesX{ nullptr };
  double *coordinatesY{ nullptr };
  yac_grid_cell yacGridCell;
};

struct CellsearchParams
{
  size_t dims[2] = { 0 };
  double (*m_xyzCoords)[3]{ nullptr };
  bool fast{ false };
};

#endif
