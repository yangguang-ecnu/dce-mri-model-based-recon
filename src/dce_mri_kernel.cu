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

/*
 * voxelConvolve
 *
 * This function represents the meat of the convolution work to be done. It is
 * provided as both a host and device function, used by the CUDA kernel but
 * also by a sequential C version for testing and benchmarking
 */
__host__ __device__ void voxelConvolve(
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
    int idx) {

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
 * gpuConvolve
 *
 * The CUDA kernel which sets up the thread index, checks boundary conditions,
 * and calls the actual work function.
 */
__global__ void gpuConvolve(
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

  voxelConvolve(
      KTrans,
      k_ep,
      dt_i,
      Ti,
      dt_j,
      Tj,
      Cpi,
      samplingRate,
      imgSeqR,
      imgSeqI,
      dimXY,
      idx);
}

/*
 * gpuSetupAndConvolve
 *
 * This is the function which is externally visible, allowing access to the
 * CUDA kernel. Its purpose is to setup the data and make the kernel call for the CUDA
 * computation, and then copy result data back out of the GPU.
 */
extern "C" float gpuSetupAndConvolve(
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

  cudaError_t errCode;

  /* Cuda memory initialization */
  errCode = cudaMalloc((void **) &d_KTrans, dimX * dimY * sizeof(float));
  if (errCode != cudaSuccess) { return -1; }
  errCode = cudaMalloc((void **) &d_k_ep, dimX * dimY * sizeof(float));
  if (errCode != cudaSuccess) { return -1; }
  errCode = cudaMalloc((void **) &d_Cpi, Ti * sizeof(float));
  if (errCode != cudaSuccess) { return -1; }
  errCode = cudaMalloc((void **) &d_imgSeqR, dimX * dimY * Tj * sizeof(float));
  if (errCode != cudaSuccess) { return -1; }
  errCode = cudaMalloc((void **) &d_imgSeqI, dimX * dimY * Tj * sizeof(float));
  if (errCode != cudaSuccess) { return -1; }

  errCode = cudaMemcpy(d_KTrans, KTrans, dimX * dimY * sizeof(float), cudaMemcpyHostToDevice);
  if (errCode != cudaSuccess) { return -1; }
  errCode = cudaMemcpy(d_k_ep, k_ep, dimX * dimY * sizeof(float), cudaMemcpyHostToDevice);
  if (errCode != cudaSuccess) { return -1; }
  errCode = cudaMemcpy(d_Cpi, Cpi, Ti * sizeof(float), cudaMemcpyHostToDevice);
  if (errCode != cudaSuccess) { return -1; }

  /* Start the timer */
  // cudaEventRecord(start_event, 0);

  /* Setup thread/block breakdown, we want dimX * dimY threads */
  int dimXY = dimX * dimY;
  int numBlocks = CEIL_DIV(dimXY, MAX_THREADS_PER_BLOCK);
  int threadsPerBlock = (numBlocks == 1 ? dimXY : MAX_THREADS_PER_BLOCK);
  dim3 dimGrid(numBlocks, 1, 1);
  dim3 dimBlock(threadsPerBlock, 1, 1);

  /* Call the kernel */
  gpuConvolve<<<dimGrid, dimBlock>>>(
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
  errCode = cudaMemcpy(imgSeqR, d_imgSeqR, dimX * dimY * Tj * sizeof(float), cudaMemcpyDeviceToHost);
  if (errCode != cudaSuccess) { return -1; }
  errCode = cudaMemcpy(imgSeqI, d_imgSeqI, dimX * dimY * Tj * sizeof(float), cudaMemcpyDeviceToHost);
  if (errCode != cudaSuccess) { return -1; }

  //return cuda_elapsed_time;
  return 0;
}

/*
 * cSetupAndConvolve
 *
 * This is a sequential C version of the same computation, provided for
 * performance comparisons and benchmarks.
 */
extern "C" float cSetupAndConvolve(
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

  int i;
  for (i = 0; i < dimX * dimY; i++) {
    voxelConvolve(
        KTrans, 
        k_ep, 
        dt_i,
        Ti,
        dt_j,
        Tj,
        Cpi, 
        samplingRate,
        imgSeqR, 
        imgSeqI, 
        dimX * dimY,
        i);
  }

  return 0;
}
