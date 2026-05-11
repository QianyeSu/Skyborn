/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef REMAP_H
#define REMAP_H

#include "varray.h"
#include "remap_vars.h"
#include "remap_grid.h"
#include "grid_cellsearch.h"
#include "grid_pointsearch.h"
#include "point.h"

constexpr double PI = std::numbers::pi;
constexpr double PI2 = (2.0 * PI);
constexpr double PIH = (0.5 * PI);

constexpr double TINY = 1.e-14;

enum struct RemapOption
{
  WriteWeights,
  GenerateWeights,
};

// clang-format off
struct  // RemapSearch
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
RemapSearch
// clang-format on
{
  RemapGrid *srcGrid{ nullptr };
  RemapGrid *tgtGrid{ nullptr };

  GridPointsearch gps{};
  GridCellsearch gcs{};
};

// clang-format off
struct  // RemapType
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
RemapType
// clang-format on
{
  int nused{};
  int gridID{ -1 };
  size_t gridsize{};
  size_t numMissVals{};
  RemapGrid srcGrid{};
  RemapGrid tgtGrid{};
  RemapVars vars{};
  RemapSearch search{};
};

int remap_check_mask_indices(const size_t (&indices)[4], Vmask const &mask);

void remap_set_option(RemapOption remapOption, int value);

void remap_init_grids(RemapMethod mapType, bool doExtrapolate, int gridID1, RemapGrid &srcGrid, int gridID2, RemapGrid &tgtGrid);

void remap_grid_free(RemapGrid &grid, bool removeMask = true);
void remap_grid_alloc(RemapMethod mapType, RemapGrid &grid);
void remap_search_init(RemapMethod mapType, RemapSearch &search, RemapGrid &srcGrid, RemapGrid &tgtGrid);

void remap_search_points(RemapSearch &rsearch, PointLonLat const &pointLL, KnnData &knnData);
int remap_search_square(RemapSearch &rsearch, PointLonLat const &pointLL, SquareCorners &squareCorners);

void remap_bilinear_weights(RemapSearch &rsearch, RemapVars &rv);
void remap_bicubic_weights(RemapSearch &rsearch, RemapVars &rv);
void remap_knn_weights(KnnParams const &knnParams, RemapSearch &rsearch, RemapVars &rv);
void remap_conserv_weights(RemapSearch &rsearch, RemapVars &rv);

void remap_bilinear(RemapSearch &rsearch, Field const &field1, Field &field2);
void remap_bicubic(RemapSearch &rsearch, Field const &field1, Field &field2);
void remap_knn(KnnParams const &knnParams, RemapSearch &rsearch, Field const &field1, Field &field2);
void remap_conserv(NormOpt normOpt, RemapSearch &rsearch, Field const &field1, Field &field2);

namespace remap
{

template <typename T>
void gradients(Varray<T> const &array, RemapGrid &grid, Vmask const &mask, RemapGradients &gradients);
void gradients(Field const &field, RemapGrid &grid, RemapGradients &gradients);

void stat(int remapOrder, RemapGrid &srcGrid, RemapGrid &tgtGrid, RemapVars &rv, Field const &field1, Field const &field2);

};  // namespace remap

void remap_write_data_scrip(std::string const &weightsfile, KnnParams const &knnParams, RemapSwitches const &remapSwitches,
                            RemapGrid &srcGrid, RemapGrid &tgtGrid, RemapVars &rv);
RemapSwitches remap_read_data_scrip(std::string const &weightsfile, int gridID1, int gridID2, RemapGrid &srcGrid,
                                    RemapGrid &tgtGrid, RemapVars &rv, KnnParams &knnParams);

int grid_search_square_reg2d_NN(size_t nx, size_t ny, size_t *nbrIndices, double *nbrDistance, double plat, double plon,
                                Varray<double> const &src_center_lat, Varray<double> const &src_center_lon);

int grid_search_square_reg2d(RemapGrid *srcGrid, SquareCorners &squareCorners, double plat, double plon);

std::pair<double, double> remap_find_weights(PointLonLat const &pointLL, SquareCorners const &squareCorners);

void remap_check_area(size_t grid_size, Varray<double> const &cell_area, const char *name);

template <typename T>
void remap_set_mask(Varray<T> const &v, size_t gridsize, size_t numMissVals, double mv, Vmask &mask);

void remap_set_mask(Field const &field1, size_t gridsize, size_t numMissVals, double missval, Vmask &imask);

#endif  // REMAP_H
