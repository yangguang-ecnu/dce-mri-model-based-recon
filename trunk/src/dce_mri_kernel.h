/* 
 * dce_mri_kernel.h
 *
 * A header file containing the public facing API of the CUDA kernel.
 *
 * Authors: Ben Felsted, Simon Williams, Cheng Ye
 *    Date: Spring 2011
 */

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
    int dimY);
