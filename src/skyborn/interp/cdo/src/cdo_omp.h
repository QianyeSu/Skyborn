#ifndef CDO_OMP_H
#define CDO_OMP_H

#ifdef _OPENMP
#define OPENMP3 200805
#define OPENMP4 201307
#define OPENMP45 201511
#define OPENMP5 201905
#define OPENMP51 202011
#define OPENMP52 202111

#if _OPENMP >= OPENMP3
#define HAVE_OPENMP3 1
#endif

#if _OPENMP >= OPENMP4
#define HAVE_OPENMP4 1
#endif

#if _OPENMP >= OPENMP45
#define HAVE_OPENMP45 1
#endif

#if _OPENMP >= OPENMP5
#define HAVE_OPENMP5 1
#endif

#if _OPENMP >= OPENMP51
#define HAVE_OPENMP51 1
#endif

#if _OPENMP >= OPENMP52
#define HAVE_OPENMP52 1
#endif
#endif

int cdo_omp_get_thread_num(void);
void cdo_omp_set_num_threads(int nthreads);

constexpr unsigned cdoMinLoopSize = 999999;

#endif
