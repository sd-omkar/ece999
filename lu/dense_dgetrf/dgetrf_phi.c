#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "mkl.h"
#include "mkl_lapacke.h"

#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)

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

	/* Local arrays */
        srand(26);
	int ipiv[n];
        int size_a = n*n;
        printf("N = %d  ", n);
	double *a0 = (double *)calloc(size_a, sizeof(double));

        for (i=0; i<size_a; i++)
          a0[i] = -10 + rand() % 20;

        // Initilize outside timing to reduce overhead
        mkl_mic_set_offload_report(1);
        #pragma offload_transfer target(mic)

        double start, end, diff;

        #pragma offload target(mic:0) in(n, lda) in(a0:length(size_a) ALLOC) in(ipiv:length(n) ALLOC)
        {
        start = dsecnd();
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a0, lda, ipiv);
        end = dsecnd();
        }
       
        diff = end - start;
        printf("Time = %.3lf s\t", diff);

        // Factorization Flops vaires between 2n(ku+1)kl and 2n(kl+ku+1)kl for REAL
        // Solution Flops 2n(ku + 2kl) for REAL
        //float flops = 2 * (2 * n * kl * (kl / 2 + ku + 1) + 2 * n * (ku + 2 * kl));
        //printf("GFlops/sec = %f\n", flops/(1000000000*diff));
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
