#include "remap_method_conserv.h"

void
cdo_compute_overlap_info(size_t numCells, CellSearch &cellSearch, const GridCell &tgtGridCell)
{
  auto targetCell = tgtGridCell.yacGridCell;
  auto &overlapAreas = cellSearch.partialAreas;
  auto &overlapCells = cellSearch.overlapCells;
  auto &overlapBarycenters = cellSearch.overlapBarycenters;

  // Do the clipping and get the cell for the overlapping area
  yac_cell_clipping(numCells, cellSearch.gridCells.data(), targetCell, overlapCells.data());

  // Get the partial areas for the overlapping regions
  if (overlapBarycenters == nullptr)
  {
    for (size_t i = 0; i < numCells; ++i) { overlapAreas[i] = gridcell_area(overlapCells[i]); }
  }
  else
  {
    for (size_t i = 0; i < numCells; ++i)
    {
      if (overlapCells[i].num_corners > 1)
      {
        for (int j = 0; j < 3; ++j) { overlapBarycenters[i][j] = 0.0; }
        overlapAreas[i] = yac_grid_cell_area_info(overlapCells[i], overlapBarycenters[i], 1.0);
        /*
        YAC_ASSERT((overlapBarycenters[i][0] != 0.0) || (overlapBarycenters[i][1] != 0.0)
                       || (overlapBarycenters[i][2] != 0.0),
                   "ERROR(yac_compute_overlap_info): overlap was computed, still barycenter is sphere origin");
        */
        normalise_vector(overlapBarycenters[i]);
        if (overlapAreas[i] < 0.0)
        {
          overlapAreas[i] = -overlapAreas[i];
          overlapBarycenters[i][0] = -overlapBarycenters[i][0];
          overlapBarycenters[i][1] = -overlapBarycenters[i][1];
          overlapBarycenters[i][2] = -overlapBarycenters[i][2];
        }
      }
      else { overlapAreas[i] = 0.0; }
    }
  }

#ifdef VERBOSE
  for (size_t i = 0; i < numSearchCells; ++i) cdo_print("overlap area : %lf", cellSearch.partialAreas[i]);
#endif
}

static double
get_edge_direction(const double *ref_corner, double *corner_a, double *corner_b)
{
  double edge_norm[3];
  crossproduct_kahan(corner_a, corner_b, edge_norm);
  normalise_vector(edge_norm);
  // sine of the angle between the edge and the reference corner
  double angle = edge_norm[0] * ref_corner[0] + edge_norm[1] * ref_corner[1] + edge_norm[2] * ref_corner[2];
  // if the reference corner is directly on the edge
  // (for small angles sin(x)==x)
  if (std::fabs(angle) < yac_angle_tol) return 0.0;

  return copysign(1.0, angle);
}

void
cdo_compute_concave_overlap_info(size_t numCells, CellSearch &cellSearch, const GridCell &tgtGridCell)
{
  auto targetCell = tgtGridCell.yacGridCell;
  auto &overlapAreas = cellSearch.partialAreas;
  auto &overlapCells = cellSearch.overlapCells;
  auto &sourceCells = cellSearch.gridCells;
  auto &overlapBarycenters = cellSearch.overlapBarycenters;

  // initialise barycenter coordinates if available
  if (overlapBarycenters != nullptr)
    for (size_t i = 0; i < numCells; ++i)
      for (int j = 0; j < 3; ++j) { overlapBarycenters[i][j] = 0.0; }

  // common node point to all partial target cells
  double *baseCorner = targetCell.coordinates_xyz[0];

  // the triangulation algorithm only works for cells that only have great circle edges
  enum yac_edge_type edgeTypes[3] = { YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE };
  double coordinates_xyz[3][3] = { { -1.0, -1.0, -1.0 }, { -1.0, -1.0, -1.0 }, { -1.0, -1.0, -1.0 } };
  for (int k = 0; k < 3; ++k) { coordinates_xyz[0][k] = baseCorner[k]; }

  // data structure to hold the triangles of the target cell
  yac_grid_cell partialCell;
  partialCell.array_size = 3;
  partialCell.num_corners = 3;
  partialCell.coordinates_xyz = coordinates_xyz;
  partialCell.edge_type = edgeTypes;

  // initialise overlap areas
  for (size_t i = 0; i < numCells; ++i) { overlapAreas[i] = 0.0; }

  // for all triangles of the target cell
  // (triangles a formed by first corner of the target cells and each edge of
  //  the cell; the first and last edge of the cell already, contain the
  //  first corner, therefore we can skip them)
  for (size_t cornerIdx = 1; cornerIdx < targetCell.num_corners - 1; ++cornerIdx)
  {
    auto cornerA = targetCell.coordinates_xyz[cornerIdx];
    auto cornerB = targetCell.coordinates_xyz[(cornerIdx + 1)];

    // if the current edge has a length of zero
    if (points_are_identically(cornerA, cornerB)) continue;

    auto edgeDirection = get_edge_direction(baseCorner, cornerA, cornerB);

    for (int k = 0; k < 3; ++k) { partialCell.coordinates_xyz[1][k] = cornerA[k]; }
    for (int k = 0; k < 3; ++k) { partialCell.coordinates_xyz[2][k] = cornerB[k]; }

    // clip the current target cell triangle with all source cells
    yac_cell_clipping(numCells, sourceCells.data(), partialCell, overlapCells.data());

    // Get the partial areas for the overlapping regions as sum over the partial target cells.
    for (size_t i = 0; i < numCells; ++i)
    {
      if (overlapCells[i].num_corners == 0) continue;

      if (overlapBarycenters == nullptr)
        overlapAreas[i] += gridcell_area(overlapCells[i]) * edgeDirection;
      else
        overlapAreas[i] += yac_grid_cell_area_info(overlapCells[i], overlapBarycenters[i], edgeDirection);
    }
  }

  for (size_t i = 0; i < numCells; ++i)
  {
    if (overlapAreas[i] < 0.0)
    {
      overlapAreas[i] = -overlapAreas[i];
      if (overlapBarycenters != nullptr)
      {
        overlapBarycenters[i][0] = -overlapBarycenters[i][0];
        overlapBarycenters[i][1] = -overlapBarycenters[i][1];
        overlapBarycenters[i][2] = -overlapBarycenters[i][2];
      }
    }
  }

  if (overlapBarycenters != nullptr)
    for (size_t i = 0; i < numCells; i++)
      if ((overlapAreas[i] > 0.0)
          && ((overlapBarycenters[i][0] != 0.0) || (overlapBarycenters[i][1] != 0.0) || (overlapBarycenters[i][2] != 0.0)))
      {
        normalise_vector(overlapBarycenters[i]);
      }

#ifdef VERBOSE
  for (size_t i = 0; i < numCells; ++i) cdo_print("overlap area %zu: %lf", i, overlapAreas[i]);
#endif
}
