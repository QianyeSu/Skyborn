/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef CELLSEARCH_SPHEREPART_H
#define CELLSEARCH_SPHEREPART_H

#include "cellsearch_unstruct.h"
#include "cellsearch_utils.h"
#include "varray.h"
#include "cdo_options.h"
#include "cdo_omp.h"
#include "grid_convert.h"
extern "C"
{
#include "lib/yac/src/grid_cell.h"
#include "lib/yac/src/sphere_part.h"

  void yac_get_cell_bounding_circle_reg_quad(struct yac_grid_cell quad, struct bounding_circle *bnd_circle);
  void yac_get_cell_bounding_circle_unstruct_triangle(double a[3], double b[3], double c[3], struct bounding_circle *bnd_circle);
}

class CellsearchSpherepart : public CellsearchStrategy
{
public:
  CellsearchSpherepart(size_t numCells, size_t numCorners, Varray<double> const &cornerLons, Varray<double> const &cornerLats,
                       const CellsearchParams &params)
      : m_params(params)
  {
    create(numCells, numCorners, cornerLons, cornerLats);
  }
  ~CellsearchSpherepart()
  {
    if (m_yacCellSearch) yac_bnd_sphere_part_search_delete(m_yacCellSearch);
    if (m_bndCircles) delete[] m_bndCircles;
  }

  size_t
  do_cellsearch(bool isReg2dCell, const GridCell &gridCell, Varray<size_t> &searchIndices)
  {
    auto const &yacGridCell = gridCell.yacGridCell;
    size_t numCorners = yacGridCell.num_corners;
    bounding_circle bndCircle;
    auto xyz = yacGridCell.coordinates_xyz;

    if (numCorners == 4 && isReg2dCell)
      yac_get_cell_bounding_circle_reg_quad(yacGridCell, &bndCircle);
    else if (numCorners == 3)
      yac_get_cell_bounding_circle_unstruct_triangle(xyz[0], xyz[1], xyz[2], &bndCircle);
    else
      yac_get_cell_bounding_circle(yacGridCell, &bndCircle);

    size_t numSearchCells;
    size_t *currNeighs;
    yac_bnd_sphere_part_search_do_bnd_circle_search(m_yacCellSearch, &bndCircle, 1, &currNeighs, &numSearchCells);

    if (searchIndices.size() < numSearchCells) searchIndices.resize(numSearchCells);

    size_t k = 0;
    // for (size_t i = 0; i < numSearchCells; ++i) searchIndices[i] = currNeighs[i];
    for (size_t i = 0; i < numSearchCells; ++i)
    {
      if (yac_extents_overlap(&bndCircle, &m_bndCircles[currNeighs[i]])) searchIndices[k++] = currNeighs[i];
    }
    numSearchCells = k;
    free(currNeighs);

    return numSearchCells;
  }

private:
  bounding_circle *m_bndCircles{ nullptr };
  bnd_sphere_part_search *m_yacCellSearch{ nullptr };
  const CellsearchParams &m_params;

  void
  create(size_t numCells, size_t numCorners, Varray<double> const &cornerLons, Varray<double> const &cornerLats)
  {
    Varray<enum yac_edge_type> edgeTypes(numCorners, YAC_GREAT_CIRCLE_EDGE);
    Varray<yac_grid_cell> cells(Threading::ompNumMaxThreads);
    for (int i = 0; i < Threading::ompNumMaxThreads; ++i)
    {
      cells[i].coordinates_xyz = new double[numCorners][3];
      cells[i].edge_type = edgeTypes.data();
      cells[i].num_corners = numCorners;
      cells[i].array_size = numCorners;
    }

    m_bndCircles = new bounding_circle[numCells];

#ifdef _OPENMP
#pragma omp parallel for if (numCells > cdoMinLoopSize) default(shared)
#endif
    for (size_t i = 0; i < numCells; ++i)
    {
      auto ompthID = cdo_omp_get_thread_num();
      auto &cell = cells[ompthID];
      auto xyz = cell.coordinates_xyz;

      for (size_t k = 0; k < numCorners; ++k) gcLLtoXYZ(cornerLons[i * numCorners + k], cornerLats[i * numCorners + k], xyz[k]);

      if (numCorners == 3)
        yac_get_cell_bounding_circle_unstruct_triangle(xyz[0], xyz[1], xyz[2], &m_bndCircles[i]);
      else
        yac_get_cell_bounding_circle(cell, &m_bndCircles[i]);

      if (m_params.m_xyzCoords != nullptr)
      {
        auto offset = i * numCorners;
        for (size_t k = 0; k < numCorners; ++k)
          for (size_t l = 0; l < 3; ++l) m_params.m_xyzCoords[offset + k][l] = xyz[k][l];
      }
    }

    m_yacCellSearch = yac_bnd_sphere_part_search_new(m_bndCircles, numCells);

    for (int i = 0; i < Threading::ompNumMaxThreads; ++i) delete[] cells[i].coordinates_xyz;
  }
};

#endif
