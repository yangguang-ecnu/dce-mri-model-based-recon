#include <cuda.h>

#define SIZE 16
#define BLOCKSIZE 4

extern "C" void outer_compute(int *in_arr, int *out_arr);

__device__ int compare(int a, int b) {
  if (a == b) return 1;
  return 0;
}

__global__ void compute(int *d_in,int *d_out) {
  int i;

  d_out[threadIdx.x] = 0;
  for (i=0; i<SIZE/BLOCKSIZE; i++) {
      d_out[threadIdx.x] += compare(d_in[i*BLOCKSIZE+threadIdx.x],6);
  }
}

void outer_compute(int *h_in_array, int *h_out_array) {
  int *d_in_array, *d_out_array;

  /* allocate memory for device copies, and copy input to device */
  cudaMalloc((void **) &d_in_array,SIZE*sizeof(int));
  cudaMalloc((void **) &d_out_array,BLOCKSIZE*sizeof(int));
  cudaMemcpy(d_in_array,h_in_array,SIZE*sizeof(int),cudaMemcpyHostToDevice);

  /* compute number of appearances of 8 for subset of data in each thread! */
  compute<<<1,BLOCKSIZE,0>>>(d_in_array,d_out_array);
  cudaThreadSynchronize();

  cudaMemcpy(h_out_array,d_out_array,BLOCKSIZE*sizeof(int),cudaMemcpyDeviceToHost);
}
