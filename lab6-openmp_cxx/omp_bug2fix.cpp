/* FILE: omp_bug2.cpp
   DESCRIPTION:
      Another OpenMP program with a bug.
   AUTHOR: Blaise Barney
   UPDATED: Daniel R. Reynolds (updated to C++), 1/13/2013
   UPDATED: Abu Talha (updated to C++), 02/27/2025 */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {

  int nthreads, i, tid;
  float total = 0.0;  // Initialize total outside the parallel region

  // Spawn parallel region
# pragma omp parallel private(tid) reduction(+:total)
  {
    // Obtain thread number
    tid = omp_get_thread_num();
    
    // Only master thread does this
    if (tid == 0) {
      nthreads = omp_get_num_threads();
      printf("Number of threads = %i\n", nthreads);
    }

    // Everyone does this
    printf("Thread %i is starting...\n", tid);

#   pragma omp barrier

    // do some work
#   pragma omp for schedule(dynamic,10)
    for (i=0; i<1000000; i++)
      total += i * 1.0;

    printf("Thread %i is done! Partial total = %g\n", tid, total);
  } // End of parallel region

  printf("Final Total = %g\n", total);  // Print after all threads finish

  return 0;
}

