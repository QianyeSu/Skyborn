/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cellsearch_reg2d.h"
#include "rect_grid_search.h"

void
grid_boundbox_reg2d(size_t nx, size_t ny, Varray<double> const &cornerLonsReg2d, Varray<double> const &cornerLatsReg2d,
                    double *gridBoundBox)
{
  gridBoundBox[0] = cornerLatsReg2d[0];
  gridBoundBox[1] = cornerLatsReg2d[ny];
  if (gridBoundBox[0] > gridBoundBox[1])
    {
      gridBoundBox[0] = cornerLatsReg2d[ny];
      gridBoundBox[1] = cornerLatsReg2d[0];
    }
  gridBoundBox[2] = cornerLonsReg2d[0];
  gridBoundBox[3] = cornerLonsReg2d[nx];
}

static long
set_srch_indices(size_t numSearchCells, long nx, long imin, long imax, long jmin, long jmax, Varray<size_t> &srchIndices)
{
  size_t maxSize = numSearchCells + (jmax - jmin + 1) * (imax - imin + 1);
  if (srchIndices.size() < maxSize) srchIndices.resize(maxSize);
  for (long jm = jmin; jm <= jmax; ++jm)
    for (long im = imin; im <= imax; ++im) srchIndices[numSearchCells++] = jm * nx + im;

  return numSearchCells;
}

static void
check_lon_bounds(double offset, double lonMin, double lonMax, double &boundLon1, double &boundLon2)
{
  boundLon1 += offset;
  boundLon2 += offset;
  if (boundLon1 < lonMin && boundLon2 > lonMin) boundLon1 = lonMin;
  if (boundLon2 > lonMax && boundLon1 < lonMax) boundLon2 = lonMax;
}
/*
static void
debug_message(const char *txt, long imin, long imax, long jmin, long jmax, Varray<double> const &srcCornersLons)
{
  printf("%s:  lonMin=%g lonMax=%g  iMin=%ld iMax=%ld  jMin=%ld jMax %ld  numCells=%ld\n", txt, RAD2DEG * srcCornersLons[imin],
         RAD2DEG * srcCornersLons[imax + 1], imin, imax, jmin, jmax, (jmax - jmin + 1) * (imax - imin + 1));
}
*/
static size_t
get_search_cells_reg2d(long nx, long ny, Varray<double> const &srcCornersLats, Varray<double> const &srcCornersLons,
                       double const *tgtCellBoundBox, Varray<size_t> &srchIndices)
{
  auto debug{ false };
  size_t numSearchCells = 0;  // num cells in restricted search arrays

  auto nxp1 = nx + 1;
  auto nyp1 = ny + 1;

  auto imin = nxp1, imax = -1L, jmin = nyp1, jmax = -1L;

  auto lfound = rect_grid_search2(jmin, jmax, tgtCellBoundBox[0], tgtCellBoundBox[1], nyp1, srcCornersLats);
  if (!lfound) return 0;
  // printf("lfound, jmin, jmax %d %ld %ld\n", lfound, jmin, jmax);
  // if (jmin > 0) jmin--;
  // if (jmax < (ny-2)) jmax++;

  auto srcLonMin = srcCornersLons[0];
  auto srcLonMax = srcCornersLons[nx];

  auto boundLon1 = tgtCellBoundBox[2];
  auto boundLon2 = tgtCellBoundBox[3];
  // debug = (bound_lon1 <= 0 && bound_lon2 >= 0);

  if (boundLon1 <= srcLonMax && boundLon2 >= srcLonMin)
    {
      check_lon_bounds(0.0, srcLonMin, srcLonMax, boundLon1, boundLon2);
      lfound = rect_grid_search2(imin, imax, boundLon1, boundLon2, nxp1, srcCornersLons);
      if (lfound)
        {
          // if (debug) debug_message("1", imin, imax, jmin, jmax, srcCornersLons);
          numSearchCells = set_srch_indices(numSearchCells, nx, imin, imax, jmin, jmax, srchIndices);
        }
    }

  boundLon1 = tgtCellBoundBox[2];
  boundLon2 = tgtCellBoundBox[3];
  if (boundLon1 < srcLonMin || boundLon2 > srcLonMax)
    {
      if (boundLon1 <= srcLonMin && boundLon2 >= srcLonMin)
        {
          check_lon_bounds(2.0 * M_PI, srcLonMin, srcLonMax, boundLon1, boundLon2);
          long imin2 = nxp1, imax2 = -1;
          lfound = rect_grid_search2(imin2, imax2, boundLon1, boundLon2, nxp1, srcCornersLons);
          if (lfound)
            {
              if (imax != -1 && imin2 <= imax) imin2 = imax + 1;
              if (imax != -1 && imax2 <= imax) imax2 = imax + 1;
              if (imin2 >= 0 && imax2 < nx)
                {
                  // if (debug) debug_message("2", imin2, imax2, jmin, jmax, srcCornersLons);
                  numSearchCells = set_srch_indices(numSearchCells, nx, imin2, imax2, jmin, jmax, srchIndices);
                }
            }
        }

      boundLon1 = tgtCellBoundBox[2];
      boundLon2 = tgtCellBoundBox[3];
      if (boundLon1 <= srcLonMax && boundLon2 >= srcLonMax)
        {
          check_lon_bounds(-2.0 * M_PI, srcLonMin, srcLonMax, boundLon1, boundLon2);
          long imin3 = nxp1, imax3 = -1;
          lfound = rect_grid_search2(imin3, imax3, boundLon1, boundLon2, nxp1, srcCornersLons);
          if (lfound)
            {
              if (imin != nxp1 && imin3 >= imin) imin3 = imin - 1;
              if (imax != nxp1 && imax3 >= imin) imax3 = imin - 1;
              if (imin3 >= 0 && imin3 < nx)
                {
                  // if (debug) debug_message("3", imin3, imax3, jmin, jmax, srcCornersLons);
                  numSearchCells = set_srch_indices(numSearchCells, nx, imin3, imax3, jmin, jmax, srchIndices);
                }
            }
        }
    }

  if (debug) printf(" numSearchCells: %zu\n", numSearchCells);

  return numSearchCells;
}

static void
restrict_boundbox(double const *gridBoundBox, double *boundBox)
{
  if (boundBox[0] < gridBoundBox[0] && boundBox[1] > gridBoundBox[0]) boundBox[0] = gridBoundBox[0];
  if (boundBox[1] > gridBoundBox[1] && boundBox[0] < gridBoundBox[1]) boundBox[1] = gridBoundBox[1];

  if (boundBox[2] >= gridBoundBox[3] && (boundBox[3] - 2 * M_PI) > gridBoundBox[2])
    {
      boundBox[2] -= 2 * M_PI;
      boundBox[3] -= 2 * M_PI;
    }
  if (boundBox[3] <= gridBoundBox[2] && (boundBox[2] - 2 * M_PI) < gridBoundBox[3])
    {
      boundBox[2] += 2 * M_PI;
      boundBox[3] += 2 * M_PI;
    }
}

static void
boundbox_from_corners_reg2d(GridCell const &gridCell, double *boundBox)
{
  auto const *coordinatesX = gridCell.coordinatesX;
  auto const *coordinatesY = gridCell.coordinatesY;

  auto clat1 = coordinatesY[0];
  auto clat2 = coordinatesY[2];

  boundBox[0] = (clat2 > clat1) ? clat1 : clat2;
  boundBox[1] = (clat2 > clat1) ? clat2 : clat1;
  boundBox[2] = coordinatesX[0];
  boundBox[3] = coordinatesX[1];
}

static void
boundbox_from_corners_unstruct(GridCell const &gridCell, double *boundBox)
{
  auto const *coordinatesX = gridCell.coordinatesX;
  auto const *coordinatesY = gridCell.coordinatesY;

  auto clon = coordinatesX[0];
  auto clat = coordinatesY[0];

  boundBox[0] = clat;
  boundBox[1] = clat;
  boundBox[2] = clon;
  boundBox[3] = clon;

  auto nc = gridCell.yacGridCell.num_corners;
  for (size_t j = 1; j < nc; ++j)
    {
      clon = coordinatesX[j];
      clat = coordinatesY[j];

      if (clat < boundBox[0]) boundBox[0] = clat;
      if (clat > boundBox[1]) boundBox[1] = clat;
      if (clon < boundBox[2]) boundBox[2] = clon;
      if (clon > boundBox[3]) boundBox[3] = clon;
    }

  if (std::fabs(boundBox[3] - boundBox[2]) > M_PI)
    {
      boundBox[2] = 0;
      boundBox[3] = 2 * M_PI;
    }
}

size_t
CellsearchReg2d::do_cellsearch(bool isReg2dCell, const GridCell &gridCell, Varray<size_t> &searchIndices)
{
  double tgtCellBoundBox[4];
  if (isReg2dCell)
    boundbox_from_corners_reg2d(gridCell, tgtCellBoundBox);
  else
    boundbox_from_corners_unstruct(gridCell, tgtCellBoundBox);

  restrict_boundbox(m_gridBoundbox, tgtCellBoundBox);

  auto numSearchCells = get_search_cells_reg2d(m_nx, m_ny, m_cornerLatsReg2d, m_cornerLonsReg2d, tgtCellBoundBox, searchIndices);

  if (numSearchCells == 1 && m_nx == 1 && m_ny == 1 && is_equal(m_cornerLatsReg2d[0], m_cornerLatsReg2d[1])
      && is_equal(m_cornerLonsReg2d[0], m_cornerLonsReg2d[1]))
    numSearchCells = 0;

  return numSearchCells;
}
