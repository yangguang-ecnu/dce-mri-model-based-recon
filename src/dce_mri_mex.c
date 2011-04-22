/* NOTE: All comments in this file must be block comments for mex */

/* 
 * dce_mri_mex.c
 *
 * This file contains the matlab integration code for the CUDA kernel. Matlab
 * has its own wrapper syntax and special functions to copy data and so we
 * separate out all that logic into its own file.
 *
 * The purpose of this file is to take the data we receive from the matlab
 * function call, copy it into C data, and pass it to our kernel to handle just
 * like it was a normal C program. We then translates the results from the
 * kernel back into matlab data and setup the appropriate return values.
 *
 * Authors: Ben Felsted, Simon Williams, Cheng Ye
 *    Date: Spring 2011
 */

/* System Includes */
#include <stdio.h>
#include <stdlib.h>

/* Matlab Integration Includes */
#include "mex.h"

/* Local Includes */
#include "dce_mri_constants.h"
#include "dce_mri_kernel.h"

/* 
 * mexFunction
 *
 * This is like main for a matlab integration point. When matlab calls the
 * function which is setup to associate with this file, this will be the entry
 * point.
 *
 * int nlhs      - number of LHS arguments (return arguments)
 * mxArray *plhs - an array to fill with these return arguments
 * int nrhs      - number of RHS arguments (function arguments)
 * mxArray *prhs - an array containing these function arguments
 *
 * In our specific case, we are expecting the following input arguments:
 *
 * KTrans       - array of floats [2D, real valued]
 * k_ep         - array of floats [2D, real valued]
 * dt_i         - float
 * Ti           - int
 * dt_j         - float
 * Tj           - int
 * Cpi          - vector of floats [2D, real valued]
 * samplingRate - float
 *
 * and the following output arguments:
 *
 * imgSeq - array of floats [3D, complex valued]
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Check for the correct number of input arguments */
  if (nrhs != 8) {
    mexPrintf("Wrong number of arguments, expecting ...(kTrans, kEp, t0, samplingRate)\n");
    return;
  } 

  /* Check for the correct number of output arguments */
  if (nlhs != 1) {
    mexPrintf("Wrong number of return values, expecting 1");
    return;
  } 

  /* Extract specific mxArray input structures */
  const mxArray *mx_KTrans = prhs[0];
  const mxArray *mx_k_ep = prhs[1];
  const mxArray *mx_dt_i = prhs[2];
  const mxArray *mx_Ti = prhs[3];
  const mxArray *mx_dt_j = prhs[4];
  const mxArray *mx_Tj = prhs[5];
  const mxArray *mx_Cpi = prhs[6];
  const mxArray *mx_samplingRate = prhs[7];

  /* Extract specific input data ptrs */
  float *KTrans = (float *)mxGetPr(mx_KTrans);
  float *k_ep = (float *)mxGetPr(mx_k_ep);
  float dt_i = *(float *)mxGetPr(mx_dt_i);
  int Ti = *(int *)mxGetPr(mx_Ti);
  float dt_j = *(float *)mxGetPr(mx_dt_j);
  int Tj = *(int *)mxGetPr(mx_Tj);
  float *Cpi = (float *)mxGetPr(mx_Cpi);
  float samplingRate = *(float *)mxGetPr(mx_samplingRate);

  /* Extract other useful information from input */
  const mwSize *kTransDims = mxGetDimensions(mx_KTrans);

  /* Setup output matrix */
  mwSize ndim = 3;
  mwSize *dims = (mwSize *)mxMalloc(ndim * sizeof(mwSize));
  dims[0] = kTransDims[0]; dims[1] = kTransDims[1]; dims[2] = Tj;
  mxArray *mxImgSeq = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxCOMPLEX);

  /* Extract specific output data ptrs */
  float *imgSeqR = (float *)mxGetPr(mxImgSeq);
  float *imgSeqI = (float *)mxGetPi(mxImgSeq);

  /* Call the CUDA kernel with the translated data */
  host_compute(
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
      dims[0], 
      dims[1]);

  /* Free memory */
  mxFree(dims);

  /* Set output to computed matrix */
  plhs[0] = mxImgSeq;
}

