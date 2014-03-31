/*******************************************************************************
*   Copyright(C) 2010-2013 Intel Corporation. All Rights Reserved.
*   
*   The source code, information  and  material ("Material") contained herein is
*   owned  by Intel Corporation or its suppliers or licensors, and title to such
*   Material remains  with Intel Corporation  or its suppliers or licensors. The
*   Material  contains proprietary information  of  Intel or  its  suppliers and
*   licensors. The  Material is protected by worldwide copyright laws and treaty
*   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
*   modified, published, uploaded, posted, transmitted, distributed or disclosed
*   in any way  without Intel's  prior  express written  permission. No  license
*   under  any patent, copyright  or  other intellectual property rights  in the
*   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
*   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
*   intellectual  property  rights must  be express  and  approved  by  Intel in
*   writing.
*   
*   *Third Party trademarks are the property of their respective owners.
*   
*   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
*   this  notice or  any other notice embedded  in Materials by Intel or Intel's
*   suppliers or licensors in any way.
*
********************************************************************************
*/
/*
   LAPACKE_dgesv Example.
   ======================
 
   The program computes the solution to the system of linear
   equations with a square matrix A and multiple
   right-hand sides B, where A is the coefficient matrix:
 
     6.80  -6.05  -0.45   8.32  -9.67 
    -2.11  -3.30   2.58   2.71  -5.14 
     5.66   5.36  -2.70   4.35  -7.26 
     5.97  -4.44   0.27  -7.17   6.08 
     8.23   1.08   9.04   2.14  -6.87 

   and B is the right-hand side matrix:
 
     4.02  -1.56   9.81 
     6.19   4.00  -4.09 
    -8.22  -8.67  -4.57 
    -7.57   1.75  -8.61 
    -3.03   2.86   8.99 
 
   Description.
   ============
 
   The routine solves for X the system of linear equations A*X = B, 
   where A is an n-by-n matrix, the columns of matrix B are individual 
   right-hand sides, and the columns of X are the corresponding 
   solutions.

   The LU decomposition with partial pivoting and row interchanges is 
   used to factor A as A = P*L*U, where P is a permutation matrix, L 
   is unit lower triangular, and U is upper triangular. The factored 
   form of A is then used to solve the system of equations A*X = B.

   Example Program Results.
   ========================
 
 LAPACKE_dgesv (row-major, high-level) Example Program Results

 Solution
  -0.80  -0.39   0.96
  -0.70  -0.55   0.22
   0.59   0.84   1.90
   1.32  -0.10   5.36
   0.57   0.11   4.04

 Details of LU factorization
   8.23   1.08   9.04   2.14  -6.87
   0.83  -6.94  -7.92   6.55  -3.99
   0.69  -0.67 -14.18   7.24  -5.19
   0.73   0.75   0.02 -13.82  14.19
  -0.26   0.44  -0.59  -0.34  -3.43

 Pivot indices
      5      5      3      4      5
*/
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "mkl.h"
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, double* a, int lda );
extern void print_int_vector( char* desc, int n, int* a );

/* Parameters */
//#define N 5
#define NRHS 1
//#define LDA N
#define LDB NRHS

/* Main program */
int main(int argc, char *argv[]) {
        if (argc != 3) {
          printf("Usgae %s dimension bandwidth\n", argv[0]);
          exit(1);
        }

	/* Locals */
	int n, nrhs = NRHS, lda, ldb = LDB, info;
        n = atoi(argv[1]);
        lda = n;
        int kl = atoi(argv[2]);
        int ku = kl;
        int i;
        struct timeval begin, end;

	/* Local arrays */
        srand(26);
	int ipiv[n];
        // a needs to have n columns and (2*kl + ku + 1) rows
        // Relevant data starts from (kl*n+ku)th element
        // Rest initialized to 0
        int size_a = (2*kl + ku + 1)*n;
        printf("N = %d  k = %d\t", n, kl);
	double *a0 = (double *)calloc(size_a, sizeof(double));

        for (i=0; i<size_a; i++)
          //a0[i] = a1[i] = a2[i] = a3[i] = a4[i] = a5[i] = a6[i] = a7[i] = a8[i] = a9[i] = 0.;
          a0[i] = 0;
        for (i=(kl*n+ku); i<size_a; i++)
          //a0[i] = a1[i] = a2[i] = a3[i] = a4[i] = a5[i] = a6[i] = a7[i] = a8[i] = a9[i] = -10 + rand() % 20;
          a0[i] = -10 + rand() % 20;

        // RHS
	double b[LDB*n];
        for (i=0; i<LDB*n; i++)
          b[i] = -10 + rand() % 20;

	/* Executable statements */
	//printf( "LAPACKE_dgbsv (row-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */

        //mkl_mic_set_workdivision(MKL_TARGET_MIC, 0, 1);         // Full offload to mic
        //mkl_mic_set_max_memory(MKL_TARGET_MIC, 0, 4194304);     // Use upto 4GB

        gettimeofday(&begin, NULL);

	info = LAPACKE_dgbtrf(LAPACK_ROW_MAJOR, n, n, kl, ku, a0, lda, ipiv);
  
        gettimeofday(&end, NULL);
        double diff = (end.tv_usec + 1000000 * end.tv_sec) - (begin.tv_usec + 1000000 * begin.tv_sec);
        printf("Time = %.3lf ms\n", diff/1000);

	/* Check for the exact singularity */
	if( info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
		printf( "the solution could not be computed.\n" );
		exit( 1 );
	}

        // Factorization Flops vaires between 2n(ku+1)kl and 2n(kl+ku+1)kl for REAL
        // Solution Flops 2n(ku + 2kl) for REAL
        //float flops = 2 * (2 * n * kl * (kl + ku + 1));
        //printf("GFlops/sec = %f\n", flops/(1000000000*diff/1000000));
	/* Print solution */
	//print_matrix( "Solution", n, nrhs, b, ldb );
	/* Print details of LU factorization */
	//print_matrix( "Details of LU factorization", n, n, a, lda );
	/* Print pivot indices */
	//print_int_vector( "Pivot indices", n, ipiv );
        free(a0);
	exit( 0 );
} /* End of LAPACKE_dgesv Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, int n, int* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}
