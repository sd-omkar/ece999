
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
        if (argc != 2) {
          printf("Usgae %s dimension bandwidth\n", argv[0]);
          exit(1);
        }

	/* Locals */
	int n, nrhs = NRHS, lda, ldb = LDB, info;
        n = atoi(argv[1]);
        lda = n;
        int i;
        struct timeval begin, end;

	/* Local arrays */
        srand(26);
	int ipiv[n];
        // a needs to have n columns and (2*kl + ku + 1) rows
        // Relevant data starts from (kl*n+ku)th element
        // Rest initialized to 0
        int size_a = n*n;
        printf("N = %d  ", n);
	double *a0 = (double *)calloc(size_a, sizeof(double));

        for (i=0; i<size_a; i++)
          a0[i] = -10 + rand() % 20;

        mkl_mic_disable();                                       // Enable auto offload
        //mkl_mic_set_workdivision(MKL_TARGET_MIC, 0, 1);         // Full offload to mic
        //mkl_mic_set_max_memory(MKL_TARGET_MIC, 0, 4194304);     // Use upto 4GB
        mkl_mic_set_offload_report(1);                          // Generate offload report

        gettimeofday(&begin, NULL);

	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a0, lda, ipiv);
  
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

        // Factorization Flops n^3*2/3 for REAL
        //double flops = 2 * (2 * n * n * n)/3;
        //printf("GFlops/sec = %lf\n", flops/(1000000000*diff/1000000));
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
