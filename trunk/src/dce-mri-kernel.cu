/* 
 * dce-mri-kernel.cu
 *
 * Authors: Ben Felsted, Simon Williams, Cheng Ye
 *    Date: Spring 2011
 */

/* CUDA Includes */
#include <cuda.h>

/* Local Includes */
#include "dce-mri-kernel.h"

extern "C" void host_compute(float *kTrans, float *kEp, float *t0, float *imgSeqR, float *imgSeqI);

__global__ void compute() {
  int tx = threadIdx.x;

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
void host_compute(float *kTrans, float *kEp, float *t0, float *imgSeqR, float *imgSeqI) {

  /* STUP FOR HOST SETUP */

}
