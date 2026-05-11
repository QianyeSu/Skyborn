#ifndef REMAP_VARS_H
#define REMAP_VARS_H

#include <cstdio>  // size_t
#include <array>

#include "field.h"

struct RemapGradients
{
  Varray<double> lat;
  Varray<double> lon;
  Varray<double> latLon;

  void
  init(size_t size)
  {
    lat.resize(size);
    lon.resize(size);
    latLon.resize(size);
  }

  explicit RemapGradients(size_t size) { init(size); }
  RemapGradients() {}
};

enum struct RemapMethod
{
  UNDEF,
  BILINEAR,
  BICUBIC,
  KNN,
  CONSERV,
};

enum struct SubmapType
{
  NONE,
  LAF,
  SUM,
  AVG
};

enum struct NormOpt
{
  NONE,
  DESTAREA,
  FRACAREA
};

struct RemapSwitches
{
  RemapMethod mapType{ RemapMethod::UNDEF };
  SubmapType submapType{ SubmapType::NONE };
  int numNeighbors{ 0 };
  int remapOrder{ 0 };
};

struct RemapVars
{
  long numLinksPerValue{ -1 };
  RemapMethod mapType{ RemapMethod::UNDEF };  // identifier for remapping method
  NormOpt normOpt{ NormOpt::NONE };           // option for normalization (conserv only)
  size_t numLinks{ 0 };                       // actual number of links for remapping
  size_t maxLinks{ 0 };                       // current size of link arrays
  size_t numWeights{ 0 };                     // num of weights used in remapping
  const size_t resizeIncrement{ 1024 };       // default amount to increase array size

  Varray<size_t> linksOffset;
  Varray<size_t> linksPerValue;
  Varray<size_t> srcCellIndices;  // source grid indices for each link
  Varray<size_t> tgtCellIndices;  // target grid indices for each link
  Varray<double> weights;         // map weights for each link [maxLinks*numWeights]
};

void remap_field(Field &field2, double missval, size_t gridsize2, RemapVars const &rv, Field const &field1,
                 RemapGradients &gradients);
void remap_laf(Field &field2, double missval, size_t gridsize2, RemapVars const &rv, Field const &field1);
void remap_avg(Field &field2, double missval, size_t gridsize2, RemapVars const &rv, Field const &field1);
void remap_vars_init(RemapMethod mapType, int remapOrder, RemapVars &rv);
void remap_vars_free(RemapVars &rv);
void remap_vars_check_weights(RemapVars const &rv);

#endif /* REMAP_VARS_H */
