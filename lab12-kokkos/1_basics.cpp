/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//
// MODIFIED: Daniel R. Reynolds, SMU Mathematics
// EDITED: Abu Talha, SMU Mathematics
// 
// Goals:
// * learn how to #include Kokkos
// * learn how to initialize and finalize Kokkos runtime environment
// * learn how to allocate/free data with Kokkos
// * learn how to use Kokkos::parallel_for
// * learn how to use the Kokkos timer
// * learn how to use Kokkos::parallel_reduce
//
// ************************************************************************
//@HEADER
*/

#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

// EXERCISE: Include Kokkos_Core.hpp.
//           cmath library unnecessary after.
// #include <cmath>
#include <Kokkos_Core.hpp>

void checkSizes( int &N, int &M, int &S, int &nrepeat );

int main( int argc, char* argv[] )
{
  int N = -1;         // number of rows 2^12
  int M = -1;         // number of columns 2^10
  int S = -1;         // total size 2^22
  int nrepeat = 100;  // number of repeats of the test

  // Read command line arguments.
  for ( int i = 0; i < argc; i++ ) {
    if ( ( strcmp( argv[ i ], "-N" ) == 0 ) || ( strcmp( argv[ i ], "-Rows" ) == 0 ) ) {
      N = pow( 2, atoi( argv[ ++i ] ) );
      printf( "  User N is %d\n", N );
    }
    else if ( ( strcmp( argv[ i ], "-M" ) == 0 ) || ( strcmp( argv[ i ], "-Columns" ) == 0 ) ) {
      M = pow( 2, atof( argv[ ++i ] ) );
      printf( "  User M is %d\n", M );
    }
    else if ( ( strcmp( argv[ i ], "-S" ) == 0 ) || ( strcmp( argv[ i ], "-Size" ) == 0 ) ) {
      S = pow( 2, atof( argv[ ++i ] ) );
      printf( "  User S is %d\n", S );
    }
    else if ( strcmp( argv[ i ], "-nrepeat" ) == 0 ) {
      nrepeat = atoi( argv[ ++i ] );
    }
    else if ( ( strcmp( argv[ i ], "-h" ) == 0 ) || ( strcmp( argv[ i ], "-help" ) == 0 ) ) {
      printf( "  y^T*A*x Options:\n" );
      printf( "  -Rows (-N) <int>:      exponent num, determines number of rows 2^num (default: 2^12 = 4096)\n" );
      printf( "  -Columns (-M) <int>:   exponent num, determines number of columns 2^num (default: 2^10 = 1024)\n" );
      printf( "  -Size (-S) <int>:      exponent num, determines total matrix size 2^num (default: 2^22 = 4096*1024 )\n" );
      printf( "  -nrepeat <int>:        number of repetitions (default: 100)\n" );
      printf( "  -help (-h):            print this message\n\n" );
      exit( 1 );
    }
  }

  // Check sizes.
  checkSizes( N, M, S, nrepeat );

  // EXERCISE: Initialize Kokkos runtime.
  //           Include braces to encapsulate code between initialize and finalize calls
  
  Kokkos::initialize( argc, argv );
  {
  // Allocate y, x vectors and Matrix A:
  // EXERCISE: Use Kokkos for data allocation in default memory space (replace std::malloc with Kokkos::kokkos_malloc)
  double* y = static_cast<double*>(Kokkos::kokkos_malloc(N * sizeof(double)));
  double* x = static_cast<double*>(Kokkos::kokkos_malloc(M * sizeof(double)));
  double* A = static_cast<double*>(Kokkos::kokkos_malloc(N * M * sizeof(double)));

  // Initialize y vector.
  // EXERCISE: Convert outer loop to Kokkos::parallel_for (be certain to use KOKKOS_LAMBDA instead of [=]).
  Kokkos::parallel_for("init_y", N, KOKKOS_LAMBDA(int i){
    y[ i ] = 1;
  });

  // Initialize x vector.
  // EXERCISE: Convert outer loop to Kokkos::parallel_for (be certain to use KOKKOS_LAMBDA instead of [=]).
  Kokkos::parallel_for("init_x", M, KOKKOS_LAMBDA(int i){
    x[ i ] = 1;  
  });

  // Initialize A matrix, note 2D indexing computation.
  // EXERCISE: Convert outer loop to Kokkos::parallel_for (be certain to use KOKKOS_LAMBDA instead of [=]).
  Kokkos::parallel_for("init_A", N, KOKKOS_LAMBDA(int j){
    for ( int i = 0; i < M; ++i ) {
      A[ j * M + i ] = 1;
    }
  });

  // Timer products.
  // EXERCISE: Switch over to using the Kokkos timer (remove "begin" and "end" entirely).
  Kokkos::Timer timer;

  for ( int repeat = 0; repeat < nrepeat; repeat++ ) {
    // Application: <y,Ax> = y^T*A*x
    double result = 0;

    // EXERCISE: Convert outer loop to Kokkos::parallel_reduce.
    Kokkos::parallel_reduce(N, KOKKOS_LAMBDA(int j, double &myresult){
      double temp2 = 0;
      for ( int i = 0; i < M; ++i ) {
        temp2 += A[ j * M + i ] * x[ i ];
      }
      myresult += y[ j ] * temp2;
    }, result);

    // Output result.
    if ( repeat == ( nrepeat - 1 ) ) {
      printf( "  Computed result for %d x %d is %lf\n", N, M, result );
    }

    const double solution = (double) N * (double) M;

    if ( result != solution ) {
      printf( "  Error: result( %lf ) != solution( %lf )\n", result, solution );
    }
  }

  // Calculate time.
  // EXERCISE: Extract the elapsed time from the Kokkos timer.
  double time = timer.seconds();

  // Calculate bandwidth.
  // Each matrix A row (each of length M) is read once.
  // The x vector (of length M) is read N times.
  // The y vector (of length N) is read once.
  double Gbytes = 1.0e-9 * double( sizeof(double) * ( M + M * N + N ) );

  // Print results (problem size, time and bandwidth in GB/s).
  printf( "  N( %d ) M( %d ) nrepeat ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
          N, M, nrepeat, Gbytes * 1000, time, Gbytes * nrepeat / time );

  // EXERCISE: free data from default memory space using Kokkos::kokkos_free
  Kokkos::kokkos_free(A);
  Kokkos::kokkos_free(y);
  Kokkos::kokkos_free(x);

  // EXERCISE: finalize Kokkos runtime
  }
  Kokkos::finalize();

  return 0;
}

void checkSizes( int &N, int &M, int &S, int &nrepeat ) {
  // If S is undefined and N or M is undefined, set S to 2^22 or the bigger of N and M.
  if ( S == -1 && ( N == -1 || M == -1 ) ) {
    S = pow( 2, 22 );
    if ( S < N ) S = N;
    if ( S < M ) S = M;
  }

  // If S is undefined and both N and M are defined, set S = N * M.
  if ( S == -1 ) S = N * M;

  // If both N and M are undefined, fix row length to the smaller of S and 2^10 = 1024.
  if ( N == -1 && M == -1 ) {
    if ( S > 1024 ) {
      M = 1024;
    }
    else {
      M = S;
    }
  }

  // If only M is undefined, set it.
  if ( M == -1 ) M = S / N;

  // If N is undefined, set it.
  if ( N == -1 ) N = S / M;

  printf( "  Total size S = %d N = %d M = %d\n", S, N, M );

  // Check sizes.
  if ( ( S < 0 ) || ( N < 0 ) || ( M < 0 ) || ( nrepeat < 0 ) ) {
    printf( "  Sizes must be greater than 0.\n" );
    exit( 1 );
  }

  if ( ( N * M ) != S ) {
    printf( "  N * M != S\n" );
    exit( 1 );
  }
}
