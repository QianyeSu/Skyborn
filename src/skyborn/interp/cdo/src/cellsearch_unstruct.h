/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef CELLSEARCH_H
#define CELLSEARCH_H

#include <cstdio>
#include <memory>
#include "varray.h"
#include "cellsearch_utils.h"

enum struct CellsearchUnstructMethod
{
  undefined,
  spherepart,
};

class CellsearchStrategy
{
public:
  virtual ~CellsearchStrategy() = default;
  virtual size_t do_cellsearch(bool isReg2dCell, const GridCell &gridCell, Varray<size_t> &searchIndices) = 0;
};

class CellsearchUnstruct
{
public:
  // explicit CellsearchUnstruct(std::unique_ptr<CellsearchStrategy> &&strategy = {}) : m_strategy(std::move(strategy)) {}
  CellsearchUnstruct() {}
  ~CellsearchUnstruct()
  {
    if (m_strategy) delete m_strategy;
  }
  /*
    void
    set_strategy(std::unique_ptr<PointsearchStrategy> &&strategy)
    {
      m_strategy = std::move(strategy);
    }
   */
  void
  set_strategy(CellsearchStrategy *strategy)
  {
    if (m_strategy) delete m_strategy;
    m_strategy = strategy;
  }

  size_t
  do_cellsearch(bool isReg2dCell, const GridCell &gridCell, Varray<size_t> &searchIndices)
  {
    if (m_strategy) { return m_strategy->do_cellsearch(isReg2dCell, gridCell, searchIndices); }
    std::fprintf(stderr, "CellsearchUnstruct::do_cellsearch: CellsearchStrategy not initialized!\n");
    return 0;
  }

  // private :
  // std::unique_ptr<CellsearchStrategy> m_strategy;
  CellsearchStrategy *m_strategy{ nullptr };
};

#endif
