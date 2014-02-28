/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      on for unsymmetric linear systems                               */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Department of Computer Science,                    */
/*      University of Basel, Switzerland.                               */
/*      Email: olaf.schenk@unibas.ch                                    */
/* -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mmio.h"
#include "mkl.h"

#pragma warning (disable:4996)

void sort(int *col_idx, double *a, int start, int end)
{
	int i, j, it;
	double dt;
	
	for (i=end-1; i>start; i--)
	{
		for(j=start; j<i; j++)
		{
			if (col_idx[j] > col_idx[j+1])
			{
				if (a)
				{
					dt=a[j]; 
					a[j]=a[j+1]; 
					a[j+1]=dt;
				}

				it=col_idx[j]; 
				col_idx[j]=col_idx[j+1]; 
				col_idx[j+1]=it;
			}
		}
	}
}

void coo2csr(int n, int nz, double *a, int *i_idx, int *j_idx,
	     double *csr_a, int *col_idx, int *row_start)
{
	int i, l;
	
	for (i=0; i<=n; i++) row_start[i] = 0;
	for (i=0; i<nz; i++) row_start[i_idx[i]+1]++;
	for (i=0; i<n; i++) row_start[i+1] += row_start[i];
	
	for (l=0; l<nz; l++)
	{
		i = row_start[i_idx[l]];
		csr_a[i] = a[l];
		col_idx[i] = j_idx[l];
		row_start[i_idx[l]]++;
	}
	
	for (i=n; i>0; i--) row_start[i] = row_start[i-1];
	row_start[0] = 0;

	for (i=0; i<n; i++)
		sort (col_idx, csr_a, row_start[i], row_start[i+1]);
}

static void read_mtx_and_return_csr(
	int argc, char **argv, int *mm, int *nn, int **ia, int **ja, double **aa)
{
	int ret_code;
    MM_typecode matcode;
    FILE *f;
    int m, n, nz;   
    int i, *I, *J;
    double *val;

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else    
    {
		printf("\n===================================\n");
		printf("mtx file = %s\n", argv[1]);
		if ((f = fopen(argv[1], "r")) == NULL)
		{
			printf("Could not open %s\n", argv[1]); 
            exit(1);
		}
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
	if ( mm_is_pattern(matcode) 
		|| mm_is_dense(matcode)
		|| mm_is_array(matcode) )
	{
		printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}
	

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nz)) !=0)
        exit(1);

    /* reseve memory for matrices */

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

	if (m != n) exit(1);
	
	*mm = m;	
	*nn = n;
	*ia = (int*) malloc (sizeof(int) * (n+1));
    *ja = (int*) malloc (sizeof(int) * nz);
    *aa = (double*) malloc (sizeof(double) * nz);
	coo2csr(n, nz, val, I, J, *aa, *ja, *ia);
	
	free (I);
	free (J);
	free (val);
}  

int main( int argc, char **argv ) 
{
    /* Matrix data. */
	int i, nz;
	double *x, *b, *aa;
    int n;
	int *ia, *ja;
	int  mtype = 11;        /* Real unsymmetric matrix */
	int  nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    int      maxfct, mnum, phase, error, msglvl;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    
	read_mtx_and_return_csr(argc, argv, &n, &n, &ia, &ja, &aa);
	printf("\nDone with reading mtx file.\n");
	printf("m = %d, n = %d, nz = %d\n\n", n, n, ia[n]);


/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters und initialize the solvers      */
/*     internal adress pointers. This is only necessary for the FIRST   */
/*     call of the PARDISO solver.                                      */
/* ---------------------------------------------------------------------*/
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
		iparm[i] = 0;
	}
    
	iparm[0] = 1;
	iparm[1] = 2;
	iparm[2] = mkl_get_max_threads();
        iparm[3] = 0;         /* No iterative-direct algorithm */
        iparm[4] = 0;         /* No user fill-in reducing permutation */
        iparm[5] = 0;         /* Write solution into x */
        iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[9] = 13; 
	iparm[10] = 1; 
	iparm[12] = 1;
        iparm[13] = 0;        /* Output: Number of perturbed pivots */
        iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1;       /* Output: Mflops for LU factorization */
        iparm[19] = 0;        /* Output: Numbers of CG Iterations */
        iparm[60] = 2;
		
    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
	{
        sscanf( var, "%d", &num_procs );
		iparm[2] = num_procs;
	}

	printf("Number of threads = %d\n", iparm[2]);
	    
    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */


/* -------------------------------------------------------------------- */    
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
	nz = ia[n];
	for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nz; i++) {
        ja[i] += 1;
    }

	/* RHS and solution vectors. */
    x = (double*) malloc(sizeof(double) * n);
	b = (double*) malloc(sizeof(double) * n);
	for (i = 0; i < n; i++)
	{
		b[i] = 0.0;
		x[i] = 0.0;
	}

	
/* -------------------------------------------------------------------- */    
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */ 
    phase = 11; 

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, aa, ia, ja, NULL, &nrhs,
		       iparm, &msglvl, NULL, NULL, &error);
    
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
   
/* -------------------------------------------------------------------- */    
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */    
    phase = 22;

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, aa, ia, ja, NULL, &nrhs,
		       iparm, &msglvl, NULL, NULL, &error);
   
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    phase = 33;

    iparm[7] = 2;       /* Max numbers of iterative refinement steps. */

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, aa, ia, ja, NULL, &nrhs,
		       iparm, &msglvl, b, x, &error);
   
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

    printf("\nSolve completed ...\n");
   
	
/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */ 
    phase = -1;                 /* Release internal memory. */

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, aa, ia, ja, NULL, &nrhs,
		       iparm, &msglvl, NULL, NULL, &error);

	free (x);
	free (b);
	free (ia);
	free (ja);
	free (aa);

	printf("\nMemory deallocated...\n");

    return 0;
} 

