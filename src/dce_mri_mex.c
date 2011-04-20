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
 * kTrans - array of floats [2D, real valued]
 * kEp    - array of floats [2D, real valued]
 * t0     - array of floats [2D, real valued]
 *
 * and the following output arguments:
 *
 * imgSeq - array of floats [3D, complex valued]
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Check for the correct number of input arguments */
  if (nrhs != 3) {
    mexPrintf("Wrong number of arguments, expecting ...(kTrans, kEp, t0)\n");
    return;
  } 

  /* Check for the correct number of output arguments */
  if (nlhs != 1) {
    mexPrintf("Wrong number of return values, expecting 1");
    return;
  } 

  /* Extract specific mxArray input structures */
  const mxArray *mxKTrans = prhs[0];
  const mxArray *mxKEp = prhs[1];
  const mxArray *mxT0 = prhs[2];

  /* Extract specific input data ptrs */
  float *kTrans = (float *)mxGetPr(mxKTrans);
  float *kEp = (float *)mxGetPr(mxKEp);
  float *t0 = (float *)mxGetPr(mxT0);

  /* Extract other useful information from input */
  const mwSize *kTransDims = mxGetDimensions(mxKTrans);

  /* Setup output matrix */
  mwSize ndim = 3;
  mwSize *dims = (mwSize *)mxMalloc(ndim * sizeof(mwSize));
  dims[0] = kTransDims[0]; dims[1] = kTransDims[1]; dims[2] = DIMENSION3;
  mxArray *mxImgSeq = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxCOMPLEX);

  /* Extract specific output data ptrs */
  float *imgSeqR = (float *)mxGetPr(mxImgSeq);
  float *imgSeqI = (float *)mxGetPi(mxImgSeq);

  /* Call the CUDA kernel with the translated data */
  host_compute(kTrans, kEp, t0, imgSeqR, imgSeqI, dims[0], dims[1]);

  /* Free memory */
  mxFree(dims);

  /* Set output to computed matrix */
  plhs[0] = mxImgSeq;
}

