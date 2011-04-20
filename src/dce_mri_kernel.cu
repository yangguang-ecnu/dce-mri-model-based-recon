/* 
 * dce_mri_kernel.cu
 *
 * Authors: Ben Felsted, Simon Williams, Cheng Ye
 *    Date: Spring 2011
 */

/* CUDA Includes */
#include <cuda.h>
#include <cutil.h>

/* Local Includes */
#include "dce_mri_constants.h"

extern "C" float
host_compute(
    float *kTrans, 
    float *kEp, 
    float *t0, 
    float *imgSeqR, 
    float *imgSeqI, 
    int dimX, 
    int dimY);

__global__ void compute(
    float *kTrans, 
    float *kEp, 
    float *t0, 
    float *imgSeqR,
    float *imgSeqI) {
  // int tx = threadIdx.x;

  /* STUB FOR COMPUTATION */

}

/*
 * host_compute
 *
 * This is the function which is externally visible, allowing access to the
 * CUDA kernel. Its purpose is to setup the data and make the kernel call for the CUDA
 * computation, and then copy result data back out of the GPU.
 *
 * In our particular case, kTrans, kEp, and t0 are input parameters and
 * imgSeqR, imgSeqI are output parameters corresponding to the real and
 * imaginary parts of the computed values.
 */
float host_compute(
    float *kTrans, 
    float *kEp, 
    float *t0, 
    float *imgSeqR, 
    float *imgSeqI, 
    int dimX, 
    int dimY) {

  /* Testing */
  *imgSeqR = kTrans[0] + kEp[0] + t0[0]; 
  *imgSeqI = kTrans[1] + kEp[1] + t0[1]; 

  // cudaEvent_t start_event, stop_event;
  float cuda_elapsed_time;

  float *d_kTrans;
  float *d_kEp;
  float *d_t0;
  float *d_imgSeqR;
  float *d_imgSeqI;

  // CUDA_SAFE_CALL(cudaEventCreate(&start_event));
  // CUDA_SAFE_CALL(cudaEventCreate(&stop_event));

  /* Cuda memory initialization */
  cudaMalloc((void **) &d_kTrans, dimX * dimY * sizeof(float));
  cudaMalloc((void **) &d_kEp, dimX * dimY * sizeof(float));
  cudaMalloc((void **) &d_t0, dimX * dimY * sizeof(float));
  cudaMalloc((void **) &d_imgSeqR, dimX * dimY * DIMENSION3 * sizeof(float));
  cudaMalloc((void **) &d_imgSeqI, dimX * dimY * DIMENSION3 * sizeof(float));

  cudaMemcpy(d_kTrans, kTrans, dimX * dimY * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_kEp, kEp, dimX * dimY * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_t0, t0, dimX * dimY * sizeof(float), cudaMemcpyHostToDevice);

  /* Start the timer */
  // cudaEventRecord(start_event, 0);

  /* Call the kernel */
  dim3 dimGrid(1, 1, 1);
  dim3 dimBlock(1, 1, 1);
  compute<<<dimGrid, dimBlock>>>(d_kTrans, d_kEp, d_t0, d_imgSeqR, d_imgSeqI);

  /* Stop the timer */
  // cudaEventRecord(stop_event, 0);
  // cudaEventSynchronize(stop_event);

  // CUDA_SAFE_CALL(cudaEventElapsedTime(&cuda_elapsed_time, start_event, stop_event));

  /* Copy cuda memory back to device */
  cudaMemcpy(imgSeqR, d_imgSeqR, dimX * dimY * DIMENSION3 * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(imgSeqI, d_imgSeqI, dimX * dimY * DIMENSION3 * sizeof(float), cudaMemcpyDeviceToHost);

  return cuda_elapsed_time;
}
