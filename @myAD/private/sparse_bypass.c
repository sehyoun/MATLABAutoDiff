/*****************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    sparseclean
 * Filename:    sparseclean.c
 * Programmer:  James Tursa
 * Version:     1.02
 * Date:        March 15, 2013
 * Copyright:   (c) 2013 by James Tursa, All Rights Reserved
 *
  %  This code uses the BSD License:
  %
  %  Redistribution and use in source and binary forms, with or without
  %  modification, are permitted provided that the following conditions are
  %  met:
  %
  %     * Redistributions of source code must retain the above copyright
  %       notice, this list of conditions and the following disclaimer.
  %     * Redistributions in binary form must reproduce the above copyright
  %       notice, this list of conditions and the following disclaimer in
  %       the documentation and/or other materials provided with the distribution
  %
  %  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  %  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  %  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  %  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  %  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  %  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  %  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  %  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  %  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  %  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  %  POSSIBILITY OF SUCH DAMAGE.
 *
 * SPARSECLEAN is a mex function intended to clean (remove) small values or
 * values within a range from sparse matrices. The operation can produce a newly
 * allocated matrix or operate on the variable inplace.
 *
 * Building:
 *
 * SPARSECLEAN requires that a mex routine be built (one time only). This
 * process is typically self-building the first time you call the function
 * as long as you have the files sparseclean.m and sparseclean.c in the same
 * directory somewhere on the MATLAB path. If you need to manually build
 * the mex function, here are the commands:
 *
 * >> mex -setup
 *   (then follow instructions to select a C or C++ compiler of your choice)
 * >> mex sparseclean.c
 *
 * The usage is as follows:
 *
 * Syntax
 *
 *  B = sparseclean(A [,true])
 *      Cleans a sparse matrix of 0's
 *
 *  B = sparseclean(A,tol [,true])
 *      Cleans a sparse matrix of all abs(A(i)) <= tol
 *
 *  B = sparseclean(A,nan [,true])
 *      Cleans a sparse matrix of all A(i) that are nan
 *
 *  B = sparseclean(A,nan,value [,true])
 *      Replaces all A(i) that are nan with value
 *      If value is complex, then A must also be complex
 *
 *  B = sparseclean(A,lower_tol,upper_tol [,true])
 *      Cleans a real sparse matrix of all lower_tol <= A(i) <= upper_tol
 *      Cleans a complex sparse matrix of all lower_tol <= abs(A(i)) <= upper_tol
 *
 *  Where A = A double sparse matrix
 *        tol, value, lower_tol, and upper_tol = scalar numeric values
 *        true = Forces in-place operation even if A is shared
 *
 * If B is omitted, then A is cleaned in-place. If A is shared or potentially shared
 * because it is a cell element or struct field element or class property, then an
 * error will be thrown for this in-place case unless you add the true argument.
 * That is, adding true at the end will force the in-place syntax to work regardless
 * of whether A is shared or potentially shared. But doing so risks side effects
 * of changing other variables that are sharing data memory with A!
 *
 */


/* Includes ----------------- */
#include "mex.h"
#include <omp.h>
#include <string.h>

/* Prototypes ----------------- */
mwIndex sparse_bypass(mxArray *ii, mxArray *jj, mxArray *vv, mxArray *A);

/* Gateway ----------------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwIndex m, n, nnz;
  m = (mwIndex) mxGetScalar(prhs[3]);
  n = (mwIndex) mxGetScalar(prhs[4]);
  nnz = mxGetM(prhs[2]);

  plhs[0] = mxCreateSparse(m, n, nnz, mxREAL);

  sparse_bypass((mxArray *)prhs[0], (mxArray *)prhs[1], (mxArray *)prhs[2], (mxArray *)plhs[0]);
}


mwIndex sparse_bypass(mxArray *ii, mxArray *jj, mxArray *vv, mxArray *A)
{
  double *pA;
  mwIndex *ir;
  mwIndex *jc;
  double *pii;
  double *pjj;
  double *pvv;
  mwIndex nnz;
  mwIndex iter;
  mwIndex l, counter;
  mwIndex curr_loc;

  l = mxGetN(A);
  nnz = mxGetM(vv);
  pA = mxGetPr(A);
  ir = mxGetIr(A);
  jc = mxGetJc(A);
  pii = mxGetPr(ii);
  pjj = mxGetPr(jj);
  pvv = mxGetPr(vv);

  memcpy((char *) pA, (char *) pvv, sizeof(pA[0])*nnz);

  #pragma omp parallel for simd schedule(static) if(nnz > 7000)
  for( iter=0; iter<nnz; iter++ ){
    ir[iter] = (mwIndex) pii[iter]-1;
  }

  jc[0] = 0;
  counter = 0;
  for( iter=0; iter<l; iter++ ){
    while( (pjj[counter]-1) == iter){
      counter++;
    }
    jc[iter+1] = counter;
  }

  return 0;
}
