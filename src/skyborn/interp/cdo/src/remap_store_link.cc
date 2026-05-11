/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_options.h"
#include "cdo_omp.h"
#include "remap_vars.h"
#include "remap_store_link.h"

bool
is_sorted_list(size_t n, const size_t *list)
{
  for (size_t i = 1; i < n; ++i)
    if (list[i] < list[i - 1]) return false;
  return true;
}

static int
qcompare_index(const void *a, const void *b)
{
  auto const &x = static_cast<const IndexWeight *>(a)->index;
  auto const &y = static_cast<const IndexWeight *>(b)->index;
  return (x < y) ? -1 : (x > y);
}

static int
qcompare_index4(const void *a, const void *b)
{
  auto const &x = static_cast<const IndexWeight4 *>(a)->index;
  auto const &y = static_cast<const IndexWeight4 *>(b)->index;
  return (x < y) ? -1 : (x > y);
}

static void
sort_indexWeights(size_t numWeights, IndexWeight *indexWeights)
{
  std::qsort(indexWeights, numWeights, sizeof(IndexWeight), qcompare_index);
}

static void
sort_indexWeights4(IndexWeight4 *indexWeights)
{
  std::qsort(indexWeights, 4, sizeof(IndexWeight4), qcompare_index4);
}

void
store_weightlinks(int doAlloc, size_t numWeights, const size_t *indices, const double *weights, size_t cellIndex,
                  std::vector<WeightLinks> &weightLinks)
{
  weightLinks[cellIndex].nlinks = 0;
  weightLinks[cellIndex].offset = 0;

  if (numWeights)
  {
    auto indexWeights = doAlloc ? new IndexWeight[numWeights] : weightLinks[cellIndex].indexWeights;

    for (size_t i = 0; i < numWeights; ++i)
    {
      indexWeights[i].index = indices[i];
      indexWeights[i].weight = weights[i];
    }

    if (numWeights > 1 && !is_sorted_list(numWeights, indices)) sort_indexWeights(numWeights, indexWeights);

    weightLinks[cellIndex].nlinks = numWeights;

    if (doAlloc) weightLinks[cellIndex].indexWeights = indexWeights;
  }
}

void
store_weightlinks_bicubic(const size_t *indices, double (&weights)[4][4], size_t cellIndex, std::vector<WeightLinks4> &weightLinks)
{
  weightLinks[cellIndex].nlinks = 0;
  weightLinks[cellIndex].offset = 0;

  auto indexWeights = weightLinks[cellIndex].indexWeights;

  for (int i = 0; i < 4; ++i)
  {
    indexWeights[i].index = indices[i];
    for (int k = 0; k < 4; ++k) { indexWeights[i].weight[k] = weights[i][k]; }
  }

  if (!is_sorted_list(4, indices)) sort_indexWeights4(indexWeights);

  weightLinks[cellIndex].nlinks = 4;
}

void
weight_links_to_remap_links(int doAlloc, size_t tgtGridSize, std::vector<WeightLinks> &weightLinks, RemapVars &rv)
{
  size_t nlinks = 0;
  for (size_t i = 0; i < tgtGridSize; ++i)
  {
    if (weightLinks[i].nlinks)
    {
      weightLinks[i].offset = nlinks;
      nlinks += weightLinks[i].nlinks;
    }
  }

  rv.maxLinks = nlinks;
  rv.numLinks = nlinks;

  if (nlinks)
  {
    auto numWeights = rv.numWeights;
    rv.srcCellIndices.resize(nlinks);
    rv.tgtCellIndices.resize(nlinks);
    rv.weights.resize(nlinks * numWeights);
    auto &srcCellIndices = rv.srcCellIndices;
    auto &tgtCellIndices = rv.tgtCellIndices;
    auto &weights = rv.weights;

    auto useLinksArrays = (Threading::ompNumMaxThreads > 1 && rv.numLinksPerValue == -1 && tgtGridSize > cdoMinLoopSize);
    if (useLinksArrays) rv.linksOffset.resize(tgtGridSize);
    if (useLinksArrays) rv.linksPerValue.resize(tgtGridSize);

#ifdef _OPENMP
#pragma omp parallel for if (tgtGridSize > cdoMinLoopSize) schedule(static) default(shared)
#endif
    for (size_t i = 0; i < tgtGridSize; ++i)
    {
      auto numLinks = weightLinks[i].nlinks;
      if (numLinks)
      {
        auto offset = weightLinks[i].offset;
        if (useLinksArrays) rv.linksOffset[i] = offset;
        if (useLinksArrays) rv.linksPerValue[i] = numLinks;
        IndexWeight *indexWeights = weightLinks[i].indexWeights;
        for (size_t ilink = 0; ilink < numLinks; ++ilink)
        {
          srcCellIndices[offset + ilink] = indexWeights[ilink].index;
          tgtCellIndices[offset + ilink] = i;
          weights[(offset + ilink) * numWeights] = indexWeights[ilink].weight;
        }

        if (doAlloc) delete[] weightLinks[i].indexWeights;
      }
      else
      {
        if (useLinksArrays) rv.linksOffset[i] = 0;
        if (useLinksArrays) rv.linksPerValue[i] = 0;
      }
    }

    if (!doAlloc) { delete[] weightLinks[0].indexWeights; }
  }
}

void
weight_links_4_to_remap_links(size_t gridSize, std::vector<WeightLinks4> &weightLinks, RemapVars &rv)
{
  size_t nlinks = 0;
  for (size_t i = 0; i < gridSize; ++i)
  {
    if (weightLinks[i].nlinks)
    {
      weightLinks[i].offset = nlinks;
      nlinks += weightLinks[i].nlinks;
    }
  }

  rv.maxLinks = nlinks;
  rv.numLinks = nlinks;
  if (nlinks)
  {
    rv.srcCellIndices.resize(nlinks);
    rv.tgtCellIndices.resize(nlinks);
    rv.weights.resize(4 * nlinks);
    auto &srcCellIndices = rv.srcCellIndices;
    auto &tgtCellIndices = rv.tgtCellIndices;
    auto &weights = rv.weights;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (size_t i = 0; i < gridSize; ++i)
    {
      auto numLinks = weightLinks[i].nlinks;
      if (numLinks)
      {
        auto offset = weightLinks[i].offset;
        const auto indexWeights = weightLinks[i].indexWeights;
        for (size_t ilink = 0; ilink < numLinks; ++ilink)
        {
          srcCellIndices[offset + ilink] = indexWeights[ilink].index;
          tgtCellIndices[offset + ilink] = i;
          for (size_t k = 0; k < 4; ++k) { weights[(offset + ilink) * 4 + k] = indexWeights[ilink].weight[k]; }
        }
      }
    }

    delete[] weightLinks[0].indexWeights;
  }
}

void
weight_links_alloc(size_t numNeighbors, size_t gridSize, std::vector<WeightLinks> &weightLinks)
{
  weightLinks[0].indexWeights = new IndexWeight[numNeighbors * gridSize];
  for (size_t i = 1; i < gridSize; ++i) { weightLinks[i].indexWeights = weightLinks[0].indexWeights + numNeighbors * i; }
}

void
weight_links_4_alloc(size_t gridSize, std::vector<WeightLinks4> &weightLinks)
{
  weightLinks[0].indexWeights = new IndexWeight4[4 * gridSize];
  for (size_t i = 1; i < gridSize; ++i) { weightLinks[i].indexWeights = weightLinks[0].indexWeights + 4 * i; }
}
