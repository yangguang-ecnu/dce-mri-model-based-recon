#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

#define SIZE 16
#define BLOCKSIZE 4

void outer_compute(int *in_arr, int *out_arr);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int *in_array, *out_array;
  int sum=0;
  int i;

  /* initialization */
  in_array = (int *) malloc(SIZE*sizeof(int)); 
  for (i=0; i<SIZE; i++) {
    in_array[i] = random()%10; 
    printf("in_array[%d] = %d\n",i,in_array[i]);    
  }
  out_array = (int *) malloc(BLOCKSIZE*sizeof(int)); 

  /* compute number of appearances of 6 */
  outer_compute(in_array, out_array);

  for (i=0; i<BLOCKSIZE; i++) {
    sum+=out_array[i];
  }

  printf ("The number 6 appears %d times in array of  %d numbers\n",sum,SIZE);
}

