/*! \file kdtree.h
 * \brief Include file for the kdtree library
 */

/*
  Code from: https://github.com/joergdietrich/libkdtree
  Uwe Schulzweida, 20150720: changed data pointer to unsigned int
                             changed *point to point[KD_MAX_DIM]
                             changed *location to location[KD_MAX_DIM]
                             changed *min to min[KD_MAX_DIM]
                             changed *max to max[KD_MAX_DIM]
                             _compPoints: compare index if points[axis] are equal
                   20171102: renamed kd_buildArg() to kd_initArg(), changed interface and memory handling
                             changed data pointer to size_t
                   20171112: added node pool
*/
#ifndef  KDTREE_H_
#define  KDTREE_H_

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <atomic>

//#define  MEMPOOL   1

#define  KDDATATYPE     double
#define  KDATA_ABS(x)   fabs(x)

//#define  KDDATATYPE     float
//#define  KDATA_ABS(x)   fabsf(x)

typedef KDDATATYPE kdata_t;


#define KD_MAX_DIM 3


typedef struct kd_point {
    kdata_t point[KD_MAX_DIM];
    size_t index;
} kd_point;


static inline
int qcmp(struct kd_point *a, struct kd_point *b, int axis)
{
  int ret = (a->point[axis] > b->point[axis]) ? 1 : (a->point[axis] < b->point[axis]) ? -1 : 0;
  if ( ret == 0 ) ret = (a->index > b->index) ? 1 : (a->index < b->index) ? -1 : 0;

  return ret;
}

/*!
 * \struct kdNode
 * \brief kd-tree node structure definition
 */
typedef struct kdNode {
    struct kdNode *left;          /*!<the left child of the tree node */
    struct kdNode *right;         /*!<the right child of the tree node */
    kdata_t location[KD_MAX_DIM]; /*!<vector to the node's location */
    kdata_t min[KD_MAX_DIM];      /*!<vector to the min coordinates of the hyperrectangle */
    kdata_t max[KD_MAX_DIM];      /*!<vector to the max coordinates of the hyperrectangle */
    int split;                    /*!<axis along which the tree bifurcates */
    size_t index;                 /*!<optional index value */
} kdNode;


typedef struct kdNodePool {
  size_t size;
  std::atomic<size_t>  imp;
  kdNode *pool;
} kdNodePool_t;


typedef struct kdTree {
  kdNode *node;
  kdNodePool_t *nodepool;
} kdTree_t;


kdNodePool_t *kd_nodepool_new(size_t n);
void kd_nodepool_delete(kdNodePool_t *pool);


/*!
 * \struct resItem
 * \brief result items, member of a priority queue
 */
typedef struct resItem {
    struct kdNode *node;        /*!<pointer to a kdNode */
    kdata_t dist_sq;            /*!<distance squared as the priority */
} resItem;

/*!
 * \struct pqueue
 * \brief priority queue (min-max heap)
 */
typedef struct pqueue {
    struct resItem **d;       /*!<pointer to an array of result items */
    size_t size;              /*!<current length of the queue */
    size_t avail;             /*!<currently allocated queue elements */
    size_t step;              /*!<step size in which new elements are allocated */
} pqueue;

/*!
 * \struct kd_thread_data
 * \brief arguments passed to the threaded kd-Tree construction
 */
typedef struct kd_thread_data {
  kdNodePool_t *nodepool;
  struct kd_point *points;
  kdata_t min[KD_MAX_DIM];
  kdata_t max[KD_MAX_DIM];
  size_t nPoints;
  int max_threads;
  int depth;
  int dim;
} kd_thread_data;


#define KD_ORDERED   (1)
#define KD_UNORDERED (0)

/* functions for the priority queue */
struct pqueue *pqinit(struct pqueue *q, size_t n);
int pqinsert(struct pqueue *q, struct resItem *d);
struct resItem **pqremove_min(struct pqueue *q, struct resItem **d);
struct resItem **pqremove_max(struct pqueue *q, struct resItem **d);
struct resItem **pqpeek_min(struct pqueue *q, struct resItem **d);
struct resItem **pqpeek_max(struct pqueue *q, struct resItem **d);


/* general utility functions */
void *kd_malloc(size_t size, const char *msg);
int kd_isleaf(struct kdNode *n);
/* Cartesian */
kdata_t kd_sqr(kdata_t a);
kdata_t kd_min(kdata_t x, kdata_t y);
/* Spherical */
kdata_t kd_sph_orth_dist(kdata_t *p1, kdata_t *p2, int split);
kdata_t kd_sph_dist_sq(kdata_t *x, kdata_t *y);
kdata_t kd_sph_dist(kdata_t *x, kdata_t *y);
kdata_t kd_sph_bearing(kdata_t *p1, kdata_t *p2);
kdata_t kd_sph_xtd(kdata_t *p1, kdata_t *p2, kdata_t *p3);

/* helper functions for debugging */
void kd_printNode(struct kdNode *node);
void kd_printTree(struct kdNode *node);

/* Functions for building and destroying trees */
struct kdNode *kd_allocNodeP(kdNodePool *nodepool, struct kd_point *points, size_t pivot, kdata_t *min, kdata_t *max, int axis, int dim);
struct kdNode *kd_allocNode(struct kd_point *points, size_t pivot, kdata_t *min, kdata_t *max, int axis, int dim);
void kd_doDestroyTree(struct kdNode *node);
void kd_destroyTree(kdTree_t *tree);
void kd_initArg(struct kd_thread_data *d, kdNodePool_t *nodepool, struct kd_point *points, size_t nPoints,
                kdata_t *min, kdata_t *max, int depth, int max_threads, int dim);
kdTree_t *kd_buildTree(struct kd_point *points, size_t nPoints,
                            kdata_t *min, kdata_t *max, int dim, int max_threads);
void *kd_doBuildTree(void *threadarg);
kdTree_t *kd_sph_buildTree(struct kd_point *points, size_t nPoints,
                           kdata_t *min, kdata_t *max, int max_threads);
void *kd_sph_doBuildTree(void *threadarg);

/* Functions for range searches
 * Cartesian
 */
int kd_isPointInRect(struct kdNode *node, kdata_t *min, kdata_t *max, int dim);
int kd_isRectInRect(struct kdNode *node, kdata_t *min, kdata_t *max, int dim);
int kd_rectOverlapsRect(struct kdNode *node, kdata_t *min, kdata_t *max, int dim);
struct pqueue *kd_ortRangeSearch(struct kdNode *node, kdata_t *min, kdata_t *max,
                                 int dim);
int kd_doOrtRangeSearch(struct kdNode *node, kdata_t *min, kdata_t *max, int dim,
                        struct pqueue *res);
struct kdNode *kd_nearest(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq, int dim);
struct pqueue *kd_qnearest(struct kdNode *node, kdata_t *p,
                           kdata_t *max_dist_sq, size_t q, int dim);
int kd_doQnearest(struct kdNode *node, kdata_t *p,
                  kdata_t *max_dist_sq, size_t q, int dim,
                  struct pqueue *res);
struct pqueue *kd_range(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq,
                        int dim, int ordered);
int kd_doRange(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq,
               int dim, struct pqueue *res, int ordered);
/* spherical */
int kd_sph_isPointInRect(struct kdNode *node, kdata_t *min, kdata_t *max);
int kd_sph_isRectInRect(struct kdNode *node, kdata_t *min, kdata_t *max);
int kd_sph_rectOverlapsRect(struct kdNode *node, kdata_t *min, kdata_t *max);
struct pqueue *kd_sph_ortRangeSearch(struct kdNode *node, kdata_t *min,
                                     kdata_t *max);
int kd_sph_doOrtRangeSearch(struct kdNode *node, kdata_t *min, kdata_t *max,
                            struct pqueue *res);
struct kdNode *kd_sph_nearest(struct kdNode *node, kdata_t *p,
                              kdata_t *max_dist_sq);
struct pqueue *kd_sph_qnearest(struct kdNode *node, kdata_t *p,
                               kdata_t *max_dist_sq, size_t q);
int kd_sph_doQnearest(struct kdNode *node, kdata_t *p,
                      kdata_t *max_dist_sq, size_t q, struct pqueue *res);
struct pqueue *kd_sph_range(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq,
                            int ordered);
int kd_sph_doRange(struct kdNode *node, kdata_t *p, kdata_t *max_dist_sq,
                   struct pqueue *res, int ordered);

/* Functions for results heaps */
int kd_insertResTree(struct kdNode *node, struct pqueue *res);


#endif  /* KDTREE_H_ */
