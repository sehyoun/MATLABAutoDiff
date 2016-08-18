/* Inputs: prhs[0] = dA/dx corresponding to a matrix A of size (nrow x ncol)
 *         prhs[1] = nrow of matrix A
 * Outputs: plhs[0] = dA'/dx corresponding to transposed matrix A'
 * Method: This is just a row permutation of dA/dx
 *
 * by SeHyoun Ahn, Aug 2016
 */

#include <omp.h>
#include "mex.h"
#include <stdlib.h>
#include <time.h>

void insertsort(mwIndex *irs, double *prs, mwSize n) {
  mwIndex i,j;
  mwIndex swapind;
  double swapval;
  for (i=1; i<n; ++i) {
    swapind = irs[i];
    swapval = prs[i];
    for (j=i; j>=0;--j) {
      if (j==0) {
        irs[j] = swapind;
        prs[j] = swapval;
      }
      else if (swapind<irs[j-1]) {
        irs[j] = irs[j-1];
        prs[j] = prs[j-1];
      }
      else {
        irs[j] = swapind;
        prs[j] = swapval;
        break;
      }
    }
  }
};

void quicksort(mwIndex* irs, double* prs, mwSize n) {
  mwIndex front, back, pivot;
  mwIndex swapind;
  double swapval;
  pivot = rand()%n;
  front = rand()%n;
  back = rand()%n;
  if (irs[front]>irs[back]) {
    if (irs[pivot]>irs[front]) {
      pivot = irs[front];
      irs[front] = irs[0];
      irs[0] = pivot;
    }
    else if (irs[pivot]>irs[back]) {
      front = irs[pivot];
      irs[pivot] = irs[0];
      irs[0] = front;
      pivot = front;
    }
    else {
      pivot = irs[back];
      irs[back] = irs[0];
      irs[0] = pivot;
    }
  }
  else {
    if (irs[pivot]>irs[back]) {
      pivot = irs[back];
      irs[back] = irs[0];
      irs[0] = pivot;
    }
    else if (irs[pivot]>irs[front]) {
      back = irs[pivot];
      irs[pivot] = irs[0];
      irs[0] = back;
      pivot = back;
    }
    else {
      pivot = irs[front];
      irs[front] = irs[0];
      irs[0] = pivot;
    }
  }
  front = 1;
  back = n-1;
  
  while (front < back) {
    if (irs[front] < pivot) {
      ++front;
    }
    else if (irs[back] > pivot) {
      --back;
    }
    else {
      swapind = irs[back];
      swapval = prs[back];
      irs[back] = irs[front];
      prs[back] = prs[front];
      irs[front]= swapind;
      prs[front]= swapval;
      ++front;
    }
  }
  if (irs[front]<pivot) {
    swapind = irs[front];
    swapval = prs[front];
    irs[front] = irs[0];
    prs[front] = prs[0];
    irs[0] = swapind;
    prs[0] = swapval;
    if (front > 17)
      quicksort(&irs[0],&prs[0],front);
    else if (front > 1)
      insertsort(&irs[0],&prs[0],front);
    if ((n-1-front) > 17)
      quicksort(&irs[front+1],&prs[front+1],n-1-front);
    else if ((n-1-front) > 1)
      insertsort(&irs[front+1],&prs[front+1],n-1-front);
  }
  else {
    swapind = irs[front-1];
    swapval = prs[front-1];
    irs[front-1] = irs[0];
    prs[front-1] = prs[0];
    irs[0] = swapind;
    prs[0] = swapval;
    if (front-1 > 17)
      quicksort(&irs[0],&prs[0],front-1);
    else if (front-1 > 1)
      insertsort(&irs[0],&prs[0],front-1);
    if (n-front > 17)
      quicksort(&irs[front],&prs[front],n-front);
    else if (n-front > 1)
      insertsort(&irs[front],&prs[front],n-front);
  }
};

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  srand(time(NULL));
  mwSize nrow, ncol, nderiv,nnz;
  
  /* Read In Sparse Matrix */
  mwIndex *irsA, *jcsA;
  double *prA;
  nderiv = mxGetN(prhs[0]);
  ncol   = mxGetM(prhs[0]);
  irsA   = mxGetIr(prhs[0]);
  jcsA   = mxGetJc(prhs[0]);
  prA    = mxGetPr(prhs[0]);
  nnz    = jcsA[nderiv];
  
  /* Read in the number of rows in the Matrix A */
  nrow   = mxGetScalar(prhs[1]);
  ncol   = ncol/nrow;
  
  /* Prepare Output Matrix */
  mwIndex *lirs, *ljcs;
  double *lpr;
  lirs    = mxMalloc( nnz * sizeof(*lirs));
  ljcs    = mxMalloc( (nderiv+1) * sizeof(*ljcs));
  lpr     = mxMalloc( nnz * sizeof(*lpr));
  
  mwIndex i,j,tmp;
  ljcs[0]=0;
  
  #pragma omp parallel for default(shared) private(i,j) num_threads(2)
  for (i=0; i<nderiv; ++i) {
    /* Compute the new row Index */
    for (j = jcsA[i]; j<jcsA[i+1]; ++j){
      lirs[j] = (irsA[j]%nrow)*ncol + irsA[j]/nrow;
      lpr[j] = prA[j];
    }
    ljcs[i+1] = jcsA[i+1];
    /* Sort to ensure sorted CSC format */
    if ( (ljcs[i+1] - ljcs[i]) > 17)
      quicksort(&lirs[ljcs[i]], &lpr[ljcs[i]], ljcs[i+1]-ljcs[i]);
    else if ( (ljcs[i+1] - ljcs[i]) > 1)
      insertsort(&lirs[ljcs[i]], &lpr[ljcs[i]], ljcs[i+1]-ljcs[i]);
  }
  plhs[0] = mxCreateSparse(nrow*ncol,nderiv,nnz,mxREAL);
  if (nnz>0) {
    mxSetIr(plhs[0],lirs);
    mxSetJc(plhs[0],ljcs);
    mxSetPr(plhs[0],lpr);
  }
}
