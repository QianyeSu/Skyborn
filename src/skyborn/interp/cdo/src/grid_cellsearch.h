/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef GRID_CELLSEARCH_H
#define GRID_CELLSEARCH_H

#include <cstddef>
#include <string>
#include "cellsearch_reg2d.h"
#include "cellsearch_unstruct.h"
#include "cellsearch_utils.h"
#include "remap_grid.h"
#include "varray.h"

enum struct CellsearchMethod
{
  spherepart,
};

class GridCellsearch
{
public:
  GridCellsearch() {}
  ~GridCellsearch()
  {
    if (reg2d) delete reg2d;
    if (m_xyzCoords) delete[] m_xyzCoords;
  }

  CellsearchMethod method{ CellsearchMethod::spherepart };

  CellsearchParams params;

  double (*m_xyzCoords)[3]{ nullptr };

  // private:
  CellsearchReg2d *reg2d{ nullptr };
  CellsearchUnstruct unstruct;
};

void grid_cellsearch_create(GridCellsearch &gcs, const RemapGrid &remapGrid);

#endif
