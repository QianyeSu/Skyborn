/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#ifndef POINTSEARCH_UTILS_H
#define POINTSEARCH_UTILS_H

#include <cmath>
#include <algorithm>

struct PointsearchParams
{
  size_t dims[2] = { 0 };
  double searchRadius{ 0 };
  bool extrapolation{ false };
  bool isCyclic{ false };
  bool isCurve{ false };
  bool useBoundBox{ false };
};

static inline double
arc_to_chord_length(double arcLength)
{
  return 2.0 * std::sin(arcLength / 2.0);
}

static inline double
chord_to_arc_length(double chordLength)
{
  return 2.0 * std::asin(std::clamp(chordLength / 2.0, -1.0, 1.0));
}

static inline void
min_point(double *min, double *point)
{
  for (int i = 0; i < 3; ++i) min[i] = (point[i] < min[i]) ? point[i] : min[i];
}

static inline void
max_point(double *max, double *point)
{
  for (int i = 0; i < 3; ++i) max[i] = (point[i] > max[i]) ? point[i] : max[i];
}

template <typename T>
static void
adjust_bbox_min(T *min)
{
  for (int i = 0; i < 3; ++i) min[i] = (min[i] < 0) ? min[i] * 1.002 : min[i] * 0.998;
}

template <typename T>
static void
adjust_bbox_max(T *max)
{
  for (int i = 0; i < 3; ++i) max[i] = (max[i] < 0) ? max[i] * 0.998 : max[i] * 1.002;
}

#endif
