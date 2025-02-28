/* FILE: omp_getEnvInfo.cpp
   DESCRIPTION:
      OpenMP Example - Get Environment Information - C++ Version
      The master thread queries and prints selected environment information.
   AUTHOR: Blaise Barney  7/06
   UPDATED: Daniel R. Reynolds (updated to C++), 1/13/2013
   UPDATED: Abu Talha (updated to C++), 02/27/2025*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {

  // local variables
  int nthreads, tid, procs, maxt, inpar, dynamic, nested;

  // Start parallel region
# pragma omp parallel private(nthreads, tid)
  {

    // Obtain thread number
    tid = omp_get_thread_num();

    // Only master thread does this
    if (tid == 0) {
      printf("Thread %i getting environment info...\n", tid);

      // Get environment information
      procs    = omp_get_num_procs();
      nthreads = omp_get_num_threads();
      maxt     = omp_get_max_threads();
      inpar    = omp_in_parallel();
      dynamic  = omp_get_dynamic();
      nested   = omp_get_nested();

      // Print environment information
      printf("The number of processors available = %i\n", procs);
      printf("The number of threads being used = %i\n", nthreads);
      printf("The maximum number of threads available = %i\n", maxt);
      printf("If we are in a parallel region (1: Yes, 0: No) = %i\n", inpar);
      printf("If dynamic threads are enabled (1: Yes, 0: No) = %i\n", dynamic);
      printf("If nested parallelism is supported (1: Yes, 0: No) = %i\n", nested);
    }

  }  // end parallel region

  return 0;
}  // end main

