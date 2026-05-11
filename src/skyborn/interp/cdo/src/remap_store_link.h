/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef REMAP_STORE_LINK_H
#define REMAP_STORE_LINK_H

#include <vector>

// Predeclarations
struct RemapVars;

struct IndexWeight
{
  size_t index{};
  double weight{};
};

struct IndexWeight4
{
  size_t index{};
  double weight[4]{};
};

struct WeightLinks
{
  size_t nlinks{};
  size_t offset{};
  IndexWeight *indexWeights{};
};

struct WeightLinks4
{
  size_t nlinks{};
  size_t offset{};
  IndexWeight4 *indexWeights{};
};

bool is_sorted_list(size_t n, const size_t *list);

void weight_links_alloc(size_t numNeighbors, size_t gridSize, std::vector<WeightLinks> &weightLinks);
void weight_links_4_alloc(size_t gridSize, std::vector<WeightLinks4> &weightLinks);
void store_weightlinks(int doAlloc, size_t numWeights, const size_t *indices, const double *weights, size_t cellIndex,
                       std::vector<WeightLinks> &weightLinks);
void store_weightlinks_bicubic(const size_t *indices, double (&weights)[4][4], size_t cellIndex,
                               std::vector<WeightLinks4> &weightLinks);
void weight_links_to_remap_links(int doAlloc, size_t gridSize, std::vector<WeightLinks> &weightLinks, RemapVars &rv);
void weight_links_4_to_remap_links(size_t gridSize, std::vector<WeightLinks4> &weightLinks, RemapVars &rv);

#endif /* REMAP_STORE_LINK */
