/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cellsearch_spherepart.h"
#include "grid_cellsearch.h"
#include "cdo_options.h"
#include <mpim_grid.h>

CellsearchMethod cellSearchMethod(CellsearchMethod::spherepart);

static void
grid_cellsearch_create_reg2d(GridCellsearch &gcs, const RemapGrid &remapGrid)
{
  auto &params = gcs.params;
  params.dims[0] = remapGrid.dims[0];
  params.dims[1] = remapGrid.dims[1];

  gcs.reg2d = new CellsearchReg2d(remapGrid.cornerLonsReg2d, remapGrid.cornerLatsReg2d, params);
}

static void
grid_cellsearch_create_unstruct(GridCellsearch &gcs, const RemapGrid &remapGrid)
{
  auto numCells = remapGrid.size;
  auto numCorners = remapGrid.numCorners;
  auto const &cornerLons = remapGrid.cornerLons;
  auto const &cornerLats = remapGrid.cornerLats;

  gcs.method = cellSearchMethod;

  // auto method = gcs.unstructMethod;
  auto &params = gcs.params;
  if (Options::fast)
    {
      if (gcs.method == CellsearchMethod::spherepart)
        {
          gcs.m_xyzCoords = new double[numCells * numCorners][3];
          params.fast = Options::fast;
          params.m_xyzCoords = gcs.m_xyzCoords;
        }
    }

  if (gcs.method == CellsearchMethod::spherepart)
    gcs.unstruct.set_strategy(new CellsearchSpherepart(numCells, numCorners, cornerLons, cornerLats, params));
}

void
grid_cellsearch_create(GridCellsearch &gcs, const RemapGrid &remapGrid)
{
  if (remapGrid.type == RemapGridType::Reg2D)
    grid_cellsearch_create_reg2d(gcs, remapGrid);
  else
    grid_cellsearch_create_unstruct(gcs, remapGrid);
}
