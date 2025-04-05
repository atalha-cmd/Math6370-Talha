/* Daniel R. Reynolds
  Edited by Abu Talha
  SMU Mathematics
  Math 6370 */

// Inclusions
#include <iostream>
#include <cmath>
#include "mpi.h"

// Description: calculates the linear residual and its averaged 2-norm (WRMS)
int linresid(double *a, double *b, double *c, double *u, double *r,
             double *res, double &norm2, int loc_N, MPI_Comm comm) {

  // Get MPI parallelism information from comm
  int nprocs, my_id;
  if (MPI_Comm_size(comm, &nprocs) != MPI_SUCCESS) {
    std::cerr << "linresid error in MPI_Comm_size\n";
    return 1;
  }
  if (MPI_Comm_rank(comm, &my_id) != MPI_SUCCESS) {
    std::cerr << "linresid error in MPI_Comm_rank\n";
    return 1;
  }

  // Initialize send/receive values
  double s_l = u[0];
  double s_r = u[loc_N-1];
  double u_l = 0.0;
  double u_r = 0.0;

  MPI_Request rreq_l, rreq_r, sreq_l, sreq_r;
  MPI_Status status;

  // phase 1: open receive channels for neighbor values
  if (my_id != 0) { // If not the first process, receive from the left
      if(MPI_Irecv(&u_l, 1, MPI_DOUBLE, my_id-1, 100, comm, &rreq_l) != MPI_SUCCESS){
          std::cerr << "linresid error in MPI_Irecv\n";
          return 1;
      }
  }
  if (my_id != nprocs-1) { // If not the last process, receive from the right
    if(MPI_Irecv(&u_r, 1, MPI_DOUBLE, my_id+1, 101, comm, &rreq_r) != MPI_SUCCESS){
          std::cerr << "linresid error in MPI_Irecv\n";
          return 1;
      }
  }

  // phase 2: send boundary values to neighbors
  if (my_id != 0) { // If not the first process, send to the left
    if(MPI_Isend(&s_l, 1, MPI_DOUBLE, my_id-1, 101, comm, &sreq_l) != MPI_SUCCESS){
          std::cerr << "linresid error in MPI_Isend\n";
          return 1;
      }
  }
  if (my_id != nprocs-1) { // If not the last process, send to the right
    if(MPI_Isend(&s_r, 1, MPI_DOUBLE, my_id+1, 100, comm, &sreq_r) != MPI_SUCCESS){
          std::cerr << "linresid error in MPI_Isend\n";
          return 1;
      }
  }

  // Step 3: Compute interior residuals while communication is in progress
  norm2 = 0.0;
  int k;
  for (k = 1; k < loc_N - 1; k++) {
    res[k] = a[k] * u[k - 1] + b[k] * u[k] + c[k] * u[k + 1] - r[k];
    norm2 += res[k] * res[k];
  }

  // Step 4: Wait for communication to complete before using u_l and u_r
  
  if (my_id != 0){
    if (MPI_Wait(&rreq_l, &status) != MPI_SUCCESS){
      std::cerr << "error in MPI_Wait\n";
      MPI_Abort(comm, 1);
    }
  }
  if (my_id != nprocs-1){
    if (MPI_Wait(&rreq_r, &status) != MPI_SUCCESS){
      std::cerr << "error in MPI_Wait\n";
      MPI_Abort(comm, 1);
    }
  }
    if (my_id != 0){
    if (MPI_Wait(&sreq_l, &status) != MPI_SUCCESS){
      std::cerr << "error in MPI_Wait\n";
      MPI_Abort(comm, 1);
    }
  }
  if (my_id != nprocs-1){
    if (MPI_Wait(&sreq_r, &status) != MPI_SUCCESS){
      std::cerr << "error in MPI_Wait\n";
      MPI_Abort(comm, 1);
    }
  }

  // Step 5: Compute residuals at left and right boundaries using received values
  res[0] = a[0] * u_l + b[0] * u[0] + c[0] * u[1] - r[0];
  norm2 += res[0] * res[0];

  k = loc_N - 1;
  res[k] = a[k] * u[k - 1] + b[k] * u[k] + c[k] * u_r - r[k];
  norm2 += res[k] * res[k];

  // Step 6: Combine local 2-norms into global averaged 2-norm
  // (this assumes that each process has the same loc_N)
  double tmp;
  if (MPI_Allreduce(&norm2, &tmp, 1, MPI_DOUBLE, MPI_SUM, comm) != MPI_SUCCESS) {
    std::cerr << "linresid error in MPI_Allreduce\n";
    return 1;
  }
  norm2 = sqrt(tmp) / loc_N / nprocs;

  return 0;
} // end linresid
