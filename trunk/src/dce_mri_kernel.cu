/* 
 * dce_mri_kernel.cu
 *
 * Authors: Ben Felsted, Simon Williams, Cheng Ye
 *    Date: Spring 2011
 */

/* CUDA Includes */
#include <cuda.h>
#include <cutil.h>
#include <math.h>

/* Local Includes */
#include "dce_mri_constants.h"

extern "C" float
host_compute(
    float *KTrans, 
    float *k_ep, 
    float dt_i, 
    int Ti, 
    float dt_j, 
    int Tj, 
    float *Cpi, 
    float samplingRate, 
    float *imgSeqR, 
    float *imgSeqI, 
    int dimX, 
    int dimY);

__global__ void compute(
    float *KTrans, 
    float *k_ep, 
    float dt_i, 
    int Ti, 
    float dt_j, 
    int Tj, 
    float *Cpi, 
    float samplingRate, 
    float *imgSeqR, 
    float *imgSeqI, 
    int dimXY, 
    int threadsPerBlock) {

  /* Get our coordinate into the KTrans, k_ep, and t0 matrices */
  int idx = blockIdx.x * threadsPerBlock + threadIdx.x;

  /* Bail early if we are an out-of-bounds thread */
  if (idx >= dimXY) { return; }

  /* Interval length */
  float L = 1.0f / samplingRate;
    
  /* Common factors */
  float my_k_ep = k_ep[idx];
  float f = my_k_ep * L;
  float a = exp(f);
  float ai = 1.0f / a;
  float b = ai - 2.0f + a;
  float c = KTrans[idx] * samplingRate / (my_k_ep * my_k_ep);

  /* Compute the convolution */
  int i, j;

  for (j = 0; j < Tj; j++) {
    imgSeqR[(dimXY * j) + idx] = 0;
  }

  for (i = 0; i < Ti; i++) {
    /* Scale the input function (vector) for the convolution */
    float ci = c * Cpi[i];
    for (j = 0; j < Tj; j++) {
      float tj = dt_j * j;
      float ti = dt_i * i;
      float u = tj - ti;

      /* More common terms */
      float g = my_k_ep * u;
      float e = exp(-g);

      /* Fake-branch (all threads take same branch) */
      float s = 0;
      if (u <= -L) {
        s = 0;
      } else if (u <= 0) {
        s = e * ai - 1 + f + g; 
      } else if (u <= L) {
        s = e * (ai - 2) + 1 + f - g;
      } else {
        s = e * b;
      }

      /* Accumulate (update time point j) */
      /* This should be shared memory and tiled to really push speed */
      imgSeqR[(dimXY * j) + idx] += ci * s;
    }
  }
}

/*
 * host_compute
 *
 * This is the function which is externally visible, allowing access to the
 * CUDA kernel. Its purpose is to setup the data and make the kernel call for the CUDA
 * computation, and then copy result data back out of the GPU.
 *
 * In our particular case, KTrans, k_ep, ... are input parameters and
 * imgSeqR, imgSeqI are output parameters corresponding to the real and
 * imaginary parts of the computed values.
 */
float host_compute(
    float *KTrans, 
    float *k_ep, 
    float dt_i, 
    int Ti, 
    float dt_j, 
    int Tj, 
    float *Cpi, 
    float samplingRate, 
    float *imgSeqR, 
    float *imgSeqI, 
    int dimX, 
    int dimY) {

  // cudaEvent_t start_event, stop_event;
  // float cuda_elapsed_time;

  float *d_KTrans;
  float *d_k_ep;
  float *d_Cpi;
  float *d_imgSeqR;
  float *d_imgSeqI;

  // CUDA_SAFE_CALL(cudaEventCreate(&start_event));
  // CUDA_SAFE_CALL(cudaEventCreate(&stop_event));

  /* Cuda memory initialization */
  cudaMalloc((void **) &d_KTrans, dimX * dimY * sizeof(float));
  cudaMalloc((void **) &d_k_ep, dimX * dimY * sizeof(float));
  cudaMalloc((void **) &d_Cpi, Ti * sizeof(float));
  cudaMalloc((void **) &d_imgSeqR, dimX * dimY * Tj * sizeof(float));
  cudaMalloc((void **) &d_imgSeqI, dimX * dimY * Tj * sizeof(float));

  cudaMemcpy(d_KTrans, KTrans, dimX * dimY * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_k_ep, k_ep, dimX * dimY * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Cpi, Cpi, Ti * sizeof(float), cudaMemcpyHostToDevice);

  /* Start the timer */
  // cudaEventRecord(start_event, 0);

  /* Setup thread/block breakdown, we want dimX * dimY threads */
  int dimXY = dimX * dimY;
  int numBlocks = CEIL_DIV(dimXY, MAX_THREADS_PER_BLOCK);
  int threadsPerBlock = (numBlocks == 1 ? dimXY : MAX_THREADS_PER_BLOCK);
  dim3 dimGrid(numBlocks, 1, 1);
  dim3 dimBlock(threadsPerBlock, 1, 1);

  /* Call the kernel */
  compute<<<dimGrid, dimBlock>>>(
      d_KTrans, 
      d_k_ep, 
      dt_i,
      Ti,
      dt_j,
      Tj,
      d_Cpi, 
      samplingRate,
      d_imgSeqR, 
      d_imgSeqI, 
      dimX * dimY,
      threadsPerBlock);

  /* Stop the timer */
  // cudaEventRecord(stop_event, 0);
  // cudaEventSynchronize(stop_event);

  // CUDA_SAFE_CALL(cudaEventElapsedTime(&cuda_elapsed_time, start_event, stop_event));

  /* Copy cuda memory back to device */
  cudaMemcpy(imgSeqR, d_imgSeqR, dimX * dimY * Tj * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(imgSeqI, d_imgSeqI, dimX * dimY * Tj * sizeof(float), cudaMemcpyDeviceToHost);

  //return cuda_elapsed_time;
  return 0;
}
