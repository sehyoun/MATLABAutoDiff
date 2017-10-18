/*
 * Inputs: v = (n x 1) {double,indicator} valued {full,sparse} vector
 *         A = (n_a x m) double valued sparse matrix with n_a = {1, n}
 * Output: bsxfun(@times,v,A)
 *
 * by SeHyoun Ahn
 */
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    mwSize nrow;
    mwSize nderiv, nnz, nrow_A;
    mwIndex *irsA,*jcsA;
    double *prA;
    double *srA;
    mwIndex *lirs,*ljcs;
    mwSize counter = 0;

    /* Read In Value Vector */
    nrow = mxGetM(prhs[0]);

    /* Read In Derivative Matrix */
    nderiv = mxGetN(prhs[1]);
    irsA = mxGetIr(prhs[1]);
    jcsA = mxGetJc(prhs[1]);
    prA = mxGetPr(prhs[1]);
    nnz = jcsA[nderiv];
    nrow_A = mxGetM(prhs[1]);
    if (nrow_A==1) {
        nnz = nnz * nrow;
    }


    /* Allocate Data Matrix for Output*/
    lirs = mxMalloc( nnz * sizeof(*lirs));
    ljcs = mxMalloc( (nderiv+1) * sizeof(*ljcs));
    srA = mxMalloc( nnz * sizeof(*srA));

    if (nrow_A != 1) {
        /* Logical Vector */
        if (mxIsLogical(prhs[0])) {
            bool *prV;
            prV = mxGetLogicals(prhs[0]);

            /* For a full vector */
            if (!mxIsSparse(prhs[0])) {
                mwIndex k,j,i;
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;

                    for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                        k = irsA[i];
                        if (prV[k]) {
                            lirs[counter]=k;
                            srA[counter] = prA[i] * prV[k];
                            ++counter;
                        }
                    }
                }
                ljcs[nderiv]=counter;
                if (counter > 0) {
                    lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                    srA = mxRealloc(srA, counter * sizeof(*srA));
                }
            }
            /* For a Sparse Vector */
            else {
                mwIndex *irsV,*jcsV,locV,i,locA;
                mwIndex k,j,tmp;
                irsV = mxGetIr(prhs[0]);
                jcsV = mxGetJc(prhs[0]);

                /* Multiply */
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;
                    locV = 0;
                    locA = jcsA[j];

                    while ((locV < jcsV[1]) && (locA<jcsA[j+1])) {
                        if (irsV[locV] == irsA[locA]) {
                            lirs[counter] = irsA[locA];
                            srA[counter] = prA[locA] * prV[locV];
                            ++counter;
                            ++locA;
                            ++locV;
                        }
                        else {
                            if (irsV[locV] > irsA[locA]) {
                                ++locA;
                            }
                            else {
                                ++locV;
                            }
                        }
                    }
                }
                ljcs[nderiv]=counter;

                if (counter > 0) {
                    lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                    srA = mxRealloc(srA, counter * sizeof(*srA));
                }
            }
        }
        /* Real Vector */
        else {
            double *prV;
            prV  = mxGetPr(prhs[0]);

            /* For a full vector */
            if (!mxIsSparse(prhs[0])) {
                mwIndex k,j,i;
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;

                    for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                        k = irsA[i];
                        if (prV[k]!=0.0) {
                            lirs[counter]=k;
                            srA[counter] = prA[i] * prV[k];
                            ++counter;
                        }
                    }
                }
                ljcs[nderiv]=counter;

                /* Reallocate Output */
                if (counter > 0) {
                    lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                    srA = mxRealloc(srA, counter * sizeof(*srA));
                }
            }
            /* For a Sparse Vector */
            else {
                mwIndex *irsV,*jcsV,locV,locA;
                mwIndex k,j,tmp;
                irsV = mxGetIr(prhs[0]);
                jcsV = mxGetJc(prhs[0]);

                /* Multiply */
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;
                    locV = 0;
                    locA = jcsA[j];

                    while ((locV < jcsV[1]) && (locA<jcsA[j+1])) {
                        if (irsV[locV] == irsA[locA]) {
                            lirs[counter] = irsA[locA];
                            srA[counter] = prA[locA] * prV[locV];
                            ++locA;
                            ++locV;
                            ++counter;
                        }
                        else {
                            if (irsV[locV] > irsA[locA]) {
                                ++locA;
                            }
                            else {
                                ++locV;
                            }
                        }
                    }
                }
                ljcs[nderiv]=counter;

                /* Reallocate Output */
                if (counter > 0) {
                    lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                    srA = mxRealloc(srA, counter * sizeof(*srA));
                }
            }
        }
    }
    else { /* bsxfun Dimension expansion */
        /* Logical Vector */
        if (mxIsLogical(prhs[0])) {
            bool *prV;
            prV = mxGetLogicals(prhs[0]);

            /* For a full vector */
            if (!mxIsSparse(prhs[0])) {
                mwIndex k,j,i;
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;

                    if (jcsA[j] < jcsA[j+1]) {
                        for (k=0; k<nrow; ++k)
                        {
                            if (prV[k]) {
                                lirs[counter] = k;
                                srA[counter] = prA[jcsA[j]] * prV[k];
                                ++counter;
                            }
                        }
                    }
                }
                ljcs[nderiv]=counter;

                if (counter > 0) {
                    lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                    srA = mxRealloc(srA, counter * sizeof(*srA));
                }
            }
            /* For a Sparse Vector */
            else {
                mwIndex *irsV,*jcsV,locV,i,locA;
                mwIndex k,j,tmp;
                irsV = mxGetIr(prhs[0]);
                jcsV = mxGetJc(prhs[0]);

                /* Multiply */
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;
                    if (jcsA[j] < jcsA[j+1]) {
                        for (k = jcsV[0]; k < jcsV[1]; ++k) {
                            lirs[counter] = irsV[k];
                            srA[counter] = prA[jcsA[j]] * prV[k];
                            ++counter;
                        }
                    }
                }
            }
            ljcs[nderiv]=counter;
            if (counter > 0) {
                lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                srA = mxRealloc(srA, counter * sizeof(*srA));
            }
        }
        /* Real Vector */
        else {
            double *prV;
            prV  = mxGetPr(prhs[0]);

            /* For a full vector */
            if (!mxIsSparse(prhs[0])) {
                mwIndex k,j,i;
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;

                    if (jcsA[j] < jcsA[j+1]) {
                        for (k = 0; k<nrow; ++k) {
                            if (prV[k]!=0.0) {
                                lirs[counter] = k;
                                srA[counter] = prA[jcsA[j]] * prV[k];
                                ++counter;
                            }
                        }
                    }
                }
                ljcs[nderiv]=counter;

                /* Reallocate Output */
                if (counter > 0) {
                    lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                    srA = mxRealloc(srA, counter * sizeof(*srA));
                }
            }
            /* For a Sparse Vector */
            else {
                mwIndex *irsV,*jcsV,locV,locA;
                mwIndex k,j,tmp;
                irsV = mxGetIr(prhs[0]);
                jcsV = mxGetJc(prhs[0]);

                /* Multiply */
                for (j=0 ; j < nderiv ; ++j) {
                    ljcs[j] = counter;

                    if (jcsA[j] < jcsA[j+1]) {
                        for (k = jcsV[0]; k < jcsV[1]; ++k) {
                            lirs[counter] = irsV[k];
                            srA[counter] = prA[jcsA[j]] * prV[k];
                            ++counter;
                        }
                    }
                }
                ljcs[nderiv]=counter;

                /* Reallocate Output */
                if (counter > 0) {
                    lirs = mxRealloc(lirs, counter * sizeof(*lirs));
                    srA = mxRealloc(srA, counter * sizeof(*srA));
                }
            }
        }
    }

    /* Set Output */
    plhs[0] = mxCreateSparse(nrow,nderiv,counter,mxREAL);
    if (counter > 0) {
        /* ugly fix for now. Will be fixed later */
        mwIndex *tmp1;
        double *aux1;
        tmp1 = mxGetIr(plhs[0]);
        mxFree(tmp1);
        tmp1 = mxGetJc(plhs[0]);
        mxFree(tmp1);
        aux1 = mxGetPr(plhs[0]);
        mxFree(aux1);

        mxSetIr(plhs[0],lirs);
        mxSetJc(plhs[0],ljcs);
        mxSetPr(plhs[0],srA);
    }
    else {
        mxFree(lirs);
        mxFree(ljcs);
        mxFree(srA);
    }
}
