/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef CELLSEARCH_REG2D_H
#define CELLSEARCH_REG2D_H

#include "cellsearch_utils.h"
#include "varray.h"

void grid_boundbox_reg2d(size_t nx, size_t ny, Varray<double> const &cornerLonsReg2d, Varray<double> const &cornerLatsReg2d,
                         double *gridBoundBox);

class CellsearchReg2d
{
public:
  CellsearchReg2d(Varray<double> const &cornerLonsReg2d, Varray<double> const &cornerLatsReg2d, const CellsearchParams &params)
  {
    m_nx = params.dims[0];
    m_ny = params.dims[1];
    create(cornerLonsReg2d, cornerLatsReg2d);
  }
  ~CellsearchReg2d() {}

  size_t do_cellsearch(bool isReg2dCell, const GridCell &gridCell, Varray<size_t> &searchIndices);

  void
  create(Varray<double> const &cornerLonsReg2d, Varray<double> const &cornerLatsReg2d)
  {
    auto nxp1 = m_nx + 1;
    auto nyp1 = m_ny + 1;

    m_cornerLonsReg2d.resize(nxp1);
    m_cornerLatsReg2d.resize(nyp1);

    varray_copy(nxp1, cornerLonsReg2d, m_cornerLonsReg2d);
    varray_copy(nyp1, cornerLatsReg2d, m_cornerLatsReg2d);

    grid_boundbox_reg2d(m_nx, m_ny, m_cornerLonsReg2d, m_cornerLatsReg2d, m_gridBoundbox);
  }

  double m_gridBoundbox[4] = { 0 };
  Varray<double> m_cornerLonsReg2d, m_cornerLatsReg2d;

private:
  size_t m_nx{ 0 };
  size_t m_ny{ 0 };
};

#endif
