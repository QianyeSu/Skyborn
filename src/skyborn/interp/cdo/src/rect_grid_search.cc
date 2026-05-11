/*
  Narrow helper extracted from upstream CDO src/interpol.cc.
  Only the regular-grid search helpers are kept here so the Skyborn vendored
  remap subset does not need the broader CDO runtime dependencies.
*/

#include "rect_grid_search.h"

#include <cmath>
#include "compare.h"

static long
find_element(double x, long nelem, Varray<double> const &v)
{
  long ii;
  long mid = 0;
  long first = 1;
  long last = nelem;

  if (v[0] < v[nelem - 1])
  {
    if (x < v[0] || x > v[nelem - 1]) return nelem;

    for (ii = 1; ii < nelem; ++ii)
    {
      mid = (first + last) >> 1;
      if (!(x < v[mid - 1] || x > v[mid])) break;
      if (x > v[mid])
        first = mid;
      else
        last = mid;
    }
  }
  else
  {
    if (x < v[nelem - 1] || x > v[0]) return nelem;

    for (ii = 1; ii < nelem; ++ii)
    {
      mid = (first + last) >> 1;
      if (!(x < v[mid] || x > v[mid - 1])) break;
      if (x < v[mid])
        first = mid;
      else
        last = mid;
    }
  }

  if (mid > 1 && is_equal(x, v[mid - 1])) mid--;

  return mid;
}

bool
rect_grid_search(size_t &ii, size_t &jj, double x, double y, size_t nxm, size_t nym, Varray<double> const &xm,
                 Varray<double> const &ym)
{
  constexpr double rtol = 1.e-12;
  auto pointFound = false;

  jj = find_element(y, (long) nym, ym);
  if (jj >= nym && std::fabs(ym[0] - y) < rtol) jj = 1;

  if (jj < nym)
  {
    ii = find_element(x, (long) nxm, xm);
    if (ii >= nxm && std::fabs(xm[0] - x) < rtol) ii = 1;
    if (ii < nxm) pointFound = true;
  }

  return pointFound;
}

bool
rect_grid_search2(long &imin, long &imax, double xmin, double xmax, long nxm, Varray<double> const &xm)
{
  auto pointFound = false;
  imin = nxm;
  imax = -1;

  auto isAscend = (xm[0] < xm[nxm - 1]);

  auto i1 = find_element(xmin, nxm, xm);
  auto i2 = find_element(xmax, nxm, xm);

  if (i1 > 0 && i1 < nxm)
  {
    pointFound = true;

    if (isAscend)
    {
      if (i1 > 1 && xmin <= xm[i1 - 1]) i1--;
      imin = i1 - 1;
      imax = i1 - 1;
    }
    else
    {
      if (i1 < nxm - 1 && xmin <= xm[i1]) i1++;
      imin = i1 - 1;
      imax = i1 - 1;
    }
  }

  if (i2 > 0 && i2 < nxm)
  {
    pointFound = true;

    if (isAscend)
    {
      if (i2 < nxm - 1 && xmax >= xm[i2]) i2++;
      imax = i2 - 1;
      if (imin == nxm) imin = imax;
    }
    else
    {
      if (i2 > 1 && xmax >= xm[i2 - 1]) i2--;
      imin = i2 - 1;
      if (imax == -1) imax = imin;
    }
  }

  return pointFound;
}
