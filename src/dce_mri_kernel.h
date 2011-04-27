/* 
 * dce_mri_kernel.h
 *
 * A header file containing the public facing API of the kernel.
 *
 * Authors: Ben Felsted, Simon Williams, Cheng Ye
 *    Date: Spring 2011
 */

/*
 * gpuSetupAndConvolve
 *
 * This is the function which is externally visible, allowing access to the
 * CUDA kernel. Its purpose is to setup the data and make the kernel call for the CUDA
 * computation, and then copy result data back out of the GPU.
 */
float gpuSetupAndConvolve(
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

/*
 * cSetupAndConvolve
 *
 * This is a sequential C version of the same computation, provided for
 * performance comparisons and benchmarks.
 */
float cSetupAndConvolve(
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
