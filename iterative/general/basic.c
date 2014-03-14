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
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
#include "mmio.h"
#include "mkl.h"

#define size 128
#define INFO 0
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
    double *aa;
    int N;
    int *ia, *ja;

    read_mtx_and_return_csr(argc, argv, &N, &N, &ia, &ja, &aa);
    printf("\nDone with reading mtx file.\n");
    printf("m = %d, n = %d, nz = %d\n\n", N, N, ia[N]);
	
        /*---------------------------------------------------------------------------
	/* Allocate storage for the ?par parameters and the solution/rhs vectors
	/*---------------------------------------------------------------------------*/

    MKL_INT ipar[size];
    double dpar[size]; 
    double *tmp = (double *)malloc((N * (2 * N + 1) + (N * (N + 9)) / 2 + 1) * sizeof(double));
    double rhs[N];
    double computed_solution[N];

        /*---------------------------------------------------------------------------
	/* Some additional variables to use with the RCI (P)FGMRES solver
	/*---------------------------------------------------------------------------*/
    
    MKL_INT itercount;
    MKL_INT RCI_request, i, ivar;
    double dvar;
    char cvar;

	/*---------------------------------------------------------------------------
	/* Initialize variables and the right hand side through matrix-vector product
	/*---------------------------------------------------------------------------*/
  ivar = N;
  cvar = 'N';
	/*---------------------------------------------------------------------------
	/* Initialize the initial guess
	/*---------------------------------------------------------------------------*/
  for (i = 0; i < N; i++)
    {
      computed_solution[i] = 1.0;
    }
	/*---------------------------------------------------------------------------
	/* Initialize the solver
	/*---------------------------------------------------------------------------*/
  dfgmres_init (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
  if (RCI_request != 0) {
    printf("Going to FAILED\n");
    goto FAILED;
  }
	/*---------------------------------------------------------------------------
	/* Set the desired parameters:
	/* LOGICAL parameters:
	/* do residual stopping test
	/* do not request for the user defined stopping test
	/* do the check of the norm of the next generated vector automatically
	/* DOUBLE PRECISION parameters
	/* set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
	/*---------------------------------------------------------------------------*/
  ipar[7] = 0;
  ipar[8] = 1;
  ipar[9] = 0;
  ipar[11] = 1;
  dpar[0] = 1.0E-3;
	/*---------------------------------------------------------------------------
	/* Check the correctness and consistency of the newly set parameters
	/*---------------------------------------------------------------------------*/
  dfgmres_check (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar,
		 tmp);
  if (RCI_request != 0) {
    printf("Going to FAILED\n");
    goto FAILED;
  }
	/*---------------------------------------------------------------------------
	/* Print the info about the RCI FGMRES method
	/*---------------------------------------------------------------------------*/
  if (INFO == 1) {
  printf ("Some info about the current run of RCI FGMRES method:\n\n");
  if (ipar[7])
    {
      printf ("As ipar[7]=%d, the automatic test for the maximal number of ",
	      ipar[7]);
      printf ("iterations will be\nperformed\n");
    }
  else
    {
      printf ("As ipar[7]=%d, the automatic test for the maximal number of ",
	      ipar[7]);
      printf ("iterations will be\nskipped\n");
    }
  printf ("+++\n");
  if (ipar[8])
    {
      printf
	("As ipar[8]=%d, the automatic residual test will be performed\n",
	 ipar[8]);
    }
  else
    {
      printf ("As ipar[8]=%d, the automatic residual test will be skipped\n",
	      ipar[8]);
    }
  printf ("+++\n");
  if (ipar[9])
    {
      printf ("As ipar[9]=%d, the user-defined stopping test will be ",
	      ipar[9]);
      printf ("requested via\nRCI_request=2\n");
    }
  else
    {
      printf ("As ipar[9]=%d, the user-defined stopping test will not be ",
	      ipar[9]);
      printf ("requested, thus,\nRCI_request will not take the value 2\n");
    }
  printf ("+++\n");
  if (ipar[10])
    {
      printf ("As ipar[10]=%d, the Preconditioned FGMRES iterations will be ",
	      ipar[10]);
      printf
	("performed, thus,\nthe preconditioner action will be requested via");
      printf ("RCI_request=3\n");
    }
  else
    {
      printf
	("As ipar[10]=%d, the Preconditioned FGMRES iterations will not ",
	 ipar[10]);
      printf ("be performed,\nthus, RCI_request will not take the value 3\n");
    }
  printf ("+++\n");
  if (ipar[11])
    {
      printf ("As ipar[11]=%d, the automatic test for the norm of the next ",
	      ipar[11]);
      printf ("generated vector is\nnot equal to zero up to rounding and ");
      printf
	("computational errors will be performed,\nthus, RCI_request will not ");
      printf ("take the value 4\n");
    }
  else
    {
      printf ("As ipar[11]=%d, the automatic test for the norm of the next ",
	      ipar[11]);
      printf ("generated vector is\nnot equal to zero up to rounding and ");
      printf
	("computational errors will be skipped,\nthus, the user-defined test ");
      printf ("will be requested via RCI_request=4\n");
    }
  printf ("+++\n\n");
  }
	/*---------------------------------------------------------------------------
	/* Compute the solution by RCI (P)FGMRES solver without preconditioning
	/* Reverse Communication starts here
	/*---------------------------------------------------------------------------*/
ONE:dfgmres (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
	/*---------------------------------------------------------------------------
	/* If RCI_request=0, then the solution was found with the required precision
	/*---------------------------------------------------------------------------*/
  if (RCI_request == 0) {
    printf("RCI_request = %d\n", RCI_request);
    printf("Going to COMPLETE\n");
    goto COMPLETE;
  }
	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------
	/* NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
	/* therefore, in C code it is required to subtract 1 from them to get C style
	/* addresses
	/*---------------------------------------------------------------------------*/
  if (RCI_request == 1)
    {
      printf("RCI_request = %d\n", RCI_request);
      mkl_dcsrgemv (&cvar, &ivar, aa, ia, ja, &tmp[ipar[21] - 1],
		    &tmp[ipar[22] - 1]);
      printf("Going to ONE\n");
      goto ONE;
    }
	/*---------------------------------------------------------------------------
	/* If RCI_request=anything else, then dfgmres subroutine failed
	/* to compute the solution vector: computed_solution[N]
	/*---------------------------------------------------------------------------*/
  else
    {
    printf("Going to FAILED\n");
      goto FAILED;
    }
	/*---------------------------------------------------------------------------
	/* Reverse Communication ends here
	/* Get the current iteration number and the FGMRES solution (DO NOT FORGET to
	/* call dfgmres_get routine as computed_solution is still containing
	/* the initial guess!)
	/*---------------------------------------------------------------------------*/
COMPLETE:dfgmres_get (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp,
	       &itercount);
  /*
     /*---------------------------------------------------------------------------
     /* Print solution vector: computed_solution[N] and the number of iterations:
     itercount
     /*--------------------------------------------------------------------------- */
  printf (" The system has been solved \n");
  printf ("\n The following solution has been obtained: \n");
  /*
  for (i = 0; i < N; i++)
    {
      printf ("computed_solution[%d]=", i);
      printf ("%e\n", computed_solution[i]);
    }
    */
  printf ("\n Number of iterations: %d\n", itercount);
  i = 1;

	/*-------------------------------------------------------------------------*/
  /* Release internal MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable MKL Memory Manager                                   */
	/*-------------------------------------------------------------------------*/
  MKL_Free_Buffers ();
  return 0;

  /*if (itercount == expected_itercount && dvar < 1.0e-14)
    {
      printf ("\nThis example has successfully PASSED through all steps of ");
      printf ("computation!\n");
      return 0;
    }
  else
    {
      printf
	("\nThis example may have FAILED as either the number of iterations ");
      printf ("differs\nfrom the expected number of iterations %d, ",
	      expected_itercount);
      printf
	("or the computed solution\ndiffers much from the expected solution ");
      printf ("(Euclidean norm is %e), or both.\n", dvar);
      return 1;
    }*/
	/*-------------------------------------------------------------------------*/
  /* Release internal MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable MKL Memory Manager                                   */
	/*-------------------------------------------------------------------------*/
FAILED:printf
    ("\nThis example FAILED as the solver has returned the ERROR ");
  printf ("code %d", RCI_request);
  MKL_Free_Buffers ();


    free (ia);
    free (ja);
    free (aa);
    free (tmp);

    printf("\nMemory deallocated...\n");

    return 0;
} 

