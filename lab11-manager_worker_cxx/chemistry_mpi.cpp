/* Daniel R. Reynolds
  Edited by Abu Talha
  SMU Mathematics
  Math 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "mpi.h"

// Prototypes
void chem_solver(double, double*, double*, double*,
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // initialize MPI
  int ierr, numprocs, myid, tag;
  MPI_Status status;
  
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
     std::cerr << "Error in calling MPI_Init\n";
     return 1;
  }
   
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    std::cerr << "Error in calling MPI_Comm_size\n";
    return 1;
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    std::cerr << "Error in calling MPI_Comm_rank\n";
    return 1;
  }
    
  // 1. set solver input parameters
  const int maxit = 1000000;
  const double lam = 1.e-2;
  const double eps = 1.e-10;
  
  double Pbuf[4], Sbuf[3];
  bool more_work;

  // 2. input the number of intervals
  if (myid == 0) {
    // === Manager Code ===
    int n;
    std::cout << "Enter the number of intervals (0 quits):\n";
    std::cin >> n;
    if (n < 1) {
      return 1;
    }

    // 3. allocate temperature and solution arrays  
    double *T = new double[n];
    double *u = new double[n];
    double *v = new double[n];
    double *w = new double[n];
    
    
    // 4. set random temperature field, initial guesses at chemical densities
    for (int i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
    for (int i=0; i<n; i++)  u[i] = 0.35;
    for (int i=0; i<n; i++)  v[i] = 0.1;
    for (int i=0; i<n; i++)  w[i] = 0.5;

    // 5. start timer
    double stime = MPI_Wtime();


    // send initial tasks to each of the worker nodes
    int numsent = 0;
    for( int i = 1; i < numprocs && i <=n; i++){
      // fill Problem buffer
      Pbuf[0] = T[i-1];
      Pbuf[1] = u[i-1];
      Pbuf[2] = v[i-1];
      Pbuf[3] = w[i-1];
      
      ierr = MPI_Send(Pbuf, 4, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
      if (ierr != MPI_SUCCESS) {
        printf("Error in calling MPI_Send in Manager = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }
      numsent++;
    }

    // obtain the workers solutions

    for(int i = 0; i < n; i++){
      ierr = MPI_Recv(Sbuf, 3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (ierr != MPI_SUCCESS) {
        printf("Error in calling MPI_Recv in Manager = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }
      // extract the solution index
      int idx = status.MPI_TAG - 1;
      int sender = status.MPI_SOURCE;
      u[idx] = Sbuf[0];
      v[idx] = Sbuf[1];
      w[idx] = Sbuf[2];

      // send task to the awaiting worker
      if (numsent < n){
        Pbuf[0] = T[numsent];
        Pbuf[1] = u[numsent];
        Pbuf[2] = v[numsent];
        Pbuf[3] = w[numsent];
        
        ierr = MPI_Send(Pbuf, 4, MPI_DOUBLE, sender , numsent + 1, MPI_COMM_WORLD);
        if (ierr != MPI_SUCCESS) {
          printf("Error in calling MPI_Send in Manager = %i\n",ierr);
          ierr = MPI_Abort(MPI_COMM_WORLD, 1);
          return 1;
        }
        numsent++;
      } else {
        ierr = MPI_Send(Pbuf, 4, MPI_DOUBLE, sender , 0, MPI_COMM_WORLD);
        if (ierr != 0) {
          printf("Error in calling MPI_Send in Manager = %i\n",ierr);
          ierr = MPI_Abort(MPI_COMM_WORLD, 1);
          return 1;
        }
      }
    }

    // 7. stop timer
    double ftime = MPI_Wtime();
    double runtime = ftime - stime;

    // 8. output solution time
    std::cout << "     runtime = " << runtime << std::endl;

    // 9. free temperature and solution arrays
    delete[] T;
    delete[] u;
    delete[] v;
    delete[] w;

    
  } else {
      // === Worker Code ===
      double T, u, v, w;
      more_work = true;

      while(more_work){
        ierr = MPI_Recv(Pbuf, 4, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (ierr != MPI_SUCCESS) {
          std::cerr << "Error in calling MPI_Recv in Worker\n";
          ierr = MPI_Abort(MPI_COMM_WORLD,1);
          return 1;
        }
        tag = status.MPI_TAG;
        if(tag == 0){
          more_work = false;
        } else {
          T = Pbuf[0];
          u = Pbuf[1];
          v = Pbuf[2];
          w = Pbuf[3];

          int its;
          double res;
          chem_solver(T, &u, &v, &w, lam, eps, maxit, &its, &res);

          if (res >= eps) {
            std::cerr << "Error in worker " << myid << ": res = " << res << std::endl;
            ierr = MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
          }

          Sbuf[0] = u;
          Sbuf[1] = v;
          Sbuf[2] = w;

          ierr = MPI_Send(Sbuf, 3, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
          if (ierr != MPI_SUCCESS) {
            std::cerr << "Error in calling MPI_Send in Worker\n";
            ierr = MPI_Abort(MPI_COMM_WORLD,1);
            return 1;
          }
        }
      } // end while
  } // end worker


  // finalize MPI
  ierr = MPI_Finalize();

  return 0;
} // end main
