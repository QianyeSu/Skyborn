/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef INTERPOL_H
#define INTERPOL_H

#include "knndata.h"

class Field;

void interpolate(Field const &field1, Field &field2);
void intgrid_bil(Field const &field1, Field &field2);
void intgrid_1nn(Field const &field1, Field &field2);
void intgrid_knn(KnnParams const &knnParams, Field const &field1, Field &field2);

constexpr double
intlin(double x, double y1, double x1, double y2, double x2)
{
  // intlin - linear interpolation
  // Uwe Schulzweida  04/05/1995
  return (y2 * (x - x1) + y1 * (x2 - x)) / (x2 - x1);
}

#endif
