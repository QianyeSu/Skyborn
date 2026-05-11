#include "cdo_omp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int
cdo_omp_get_thread_num()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

void
cdo_omp_set_num_threads(int nthreads)
{
#ifdef _OPENMP
  if (omp_get_max_threads() != nthreads) omp_set_num_threads(nthreads);
#endif
}
