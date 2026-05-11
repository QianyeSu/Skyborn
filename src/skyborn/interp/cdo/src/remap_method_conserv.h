#ifndef REMAP_METHOD_CONSERV_H
#define REMAP_METHOD_CONSERV_H

#include "grid_cellsearch.h"

extern "C"
{
#include "lib/yac/src/clipping.h"
#include "lib/yac/src/area.h"
#include "lib/yac/src/geometry.h"
}

struct CellSearch
{
  enum yac_edge_type *edgeType = nullptr;
  size_t numCorners = 0;
  size_t maxCells = 0;
  Varray<double> partialAreas;
  Varray<yac_grid_cell> gridCells;
  Varray<yac_grid_cell> overlapCells;
  double (*overlapBarycenters)[3] = nullptr;

  void
  realloc(size_t numCells, bool allocBarycenters = false)
  {
    if (numCells > maxCells)
    {
      partialAreas.resize(numCells);
      overlapCells.resize(numCells);
      gridCells.resize(numCells);

      for (size_t i = maxCells; i < numCells; ++i)
      {
        overlapCells[i].array_size = 0;
        overlapCells[i].num_corners = 0;
        overlapCells[i].edge_type = nullptr;
        overlapCells[i].coordinates_xyz = nullptr;

        gridCells[i].array_size = numCorners;
        gridCells[i].num_corners = numCorners;
        gridCells[i].edge_type = edgeType;
        gridCells[i].coordinates_xyz = new double[numCorners][3];
      }

      if (allocBarycenters) overlapBarycenters = new double[numCells][3];

      maxCells = numCells;
    }
  }

  void
  free()
  {
    for (size_t i = 0; i < maxCells; ++i)
    {
      if (overlapCells[i].array_size > 0)
      {
        if (overlapCells[i].coordinates_xyz) { std::free(overlapCells[i].coordinates_xyz); }
        if (overlapCells[i].edge_type) { std::free(overlapCells[i].edge_type); }
      }

      delete[] gridCells[i].coordinates_xyz;
    }

    varray_free(partialAreas);
    varray_free(overlapCells);
    varray_free(gridCells);

    if (overlapBarycenters) { delete[] overlapBarycenters; }
  }
};

inline double
gridcell_area(const yac_grid_cell &cell)
{
  return (cell.num_corners > 1) ? yac_grid_cell_area(cell) : 0.0;
}

void cdo_compute_overlap_info(size_t numCells, CellSearch &cellSearch, const GridCell &tgtGridCell);
void cdo_compute_concave_overlap_info(size_t numCells, CellSearch &cellSearch, const GridCell &tgtGridCell);

#endif
