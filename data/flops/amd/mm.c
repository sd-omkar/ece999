#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


int main (int argc, char *argv[]) 
{
int     tid, nthreads, i, j, k;
int     NN = atoi(argv[1]);
double  *a = (double *)malloc(sizeof(double)*NN*NN);
double  *b = (double *)malloc(sizeof(double)*NN*NN);
double  *c = (double *)malloc(sizeof(double)*NN*NN);

/*** Spawn a parallel region explicitly scoping all variables ***/
#pragma omp parallel shared(a,b,c,nthreads) private(tid,i,j,k)
  {
  tid = omp_get_thread_num();
  if (tid == 0)
    {
    nthreads = omp_get_num_threads();
    printf("Starting matrix multiple example with %d threads\n",nthreads);
    printf("Initializing matrices...\n");
    }
  /*** Initialize matrices ***/
  #pragma omp for
  for (i=0; i<NN; i++)
    for (j=0; j<NN; j++)
      a[i*NN + j]= i+j;
  #pragma omp for
  for (i=0; i<NN; i++)
    for (j=0; j<NN; j++)
      b[i*NN + j]= i*j;
  #pragma omp for
  for (i=0; i<NN; i++)
    for (j=0; j<NN; j++)
      c[i*NN + j]= 0;

  /*** Do matrix multiply sharing iterations on outer loop ***/
  /*** Display who does which iterations for demonstration purposes ***/
  if (tid == 0)
    printf("Multiplying matrices...\n");
  #pragma omp for
  for (i=0; i<NN; i++)    
    {
    for(j=0; j<NN; j++)       
      for (k=0; k<NN; k++)
        c[i*NN + j] += a[i*NN + k] * b[k*NN + j];
    }
  }   /*** End of parallel region ***/

  free(a);
  free(b);
  free(c);
}
