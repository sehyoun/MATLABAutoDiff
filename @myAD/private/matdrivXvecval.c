/* This matrix is saved in column-major order */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{

    /* Read In Matrix A */
    mwSize mA, nderiv, nnz;
    mwIndex *irsA,*jcsA;
    double *prA;
    double *srA;
    mwIndex *lirs,*ljcs;
    mwSize ncol,counter=0,nrow;
    mwIndex *pointer, rownnz, pointed;
    double *prV;

    mA = mxGetM(prhs[0]);
    nderiv = mxGetN(prhs[0]);
    irsA = mxGetIr(prhs[0]);
    jcsA = mxGetJc(prhs[0]);
    prA = mxGetPr(prhs[0]);
    nnz = jcsA[nderiv];

    /* Allocate Data Matrix for Output*/
    lirs = mxMalloc( nnz * sizeof(*lirs));
    ljcs = mxMalloc( (nderiv+1) * sizeof(*ljcs));
    srA = mxMalloc( nnz * sizeof(*srA));

    /* Throw an error if multiplication with indicator valued vector is attempted */
    if (mxIsLogical(prhs[1])) {
        mxArray *arg;
        arg = mxCreateString("AutoDiff does not support matrix-(indicator-valued vector) multiplication. Use Subsetting instead of multiplication.");
        mexCallMATLAB(0,0,1,&arg,"error");
    }

    /* Read in Vector V */
    ncol = mxGetM(prhs[1]);
    prV = mxGetPr(prhs[1]);
    nrow = mA/ncol;
    pointer = mxMalloc( nrow * sizeof(*pointer));

    /* For a full vector */
    if (!mxIsSparse(prhs[1])) {
        mwIndex i,j,k,l,m,row_pos,isinflag;

        for (j=0 ; j < nderiv ; ++j) {
            ljcs[j] = counter;

            rownnz=0;
            for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                k = irsA[i];
                row_pos = k%nrow;
                isinflag = 0;
                for (l=0; l < rownnz; ++l) {
                    if (lirs[counter+l] == row_pos) {
                        srA[counter+l] += prA[i] * prV[k / nrow];
                        isinflag = 1;
                        break;
                    }
                    if (lirs[counter+l] > row_pos) {
                        for (m = rownnz; m >l; --m) {
                            isinflag = 1;
                            srA[counter+m] = srA[counter+m-1];
                            lirs[counter+m] = lirs[counter+m-1];
                        }
                        srA[counter+l] = prA[i] * prV[k/nrow];
                        lirs[counter+l] = row_pos;
                        ++rownnz;
                        isinflag = 1;
                        break;
                    }
                }
                if (isinflag==0) {
                    srA[counter+rownnz] = prA[i]*prV[k/nrow];
                    lirs[counter+rownnz] = row_pos;
                    ++rownnz;
                }
            }
            counter += rownnz;
        }
        ljcs[nderiv]=counter;

        if (counter > 0) {
            lirs = mxRealloc(lirs, counter * sizeof(*lirs));
            srA = mxRealloc(srA, counter * sizeof(*srA));
        }
    }
    /* For a Sparse Vector */
    else {
        mwIndex *irsV,*jcsV,*pos,old,i;
        mwIndex j,k,l,m,col_pos,row_pos,isinflag;
        irsV = mxGetIr(prhs[1]);
        jcsV = mxGetJc(prhs[1]);
        pos = mxMalloc(ncol * sizeof(*pos));

        old=0;
        for (i=0; i<jcsV[1]; ++i) {
            for (; old<irsV[i]; ++old) {
                pos[old] = -1;
            }
            pos[old] = i;
            ++old;
        }
        for (i=old; i<ncol; ++i) {
            pos[i]=-1;
        }

        /* Multiply */
        for (j=0 ; j < nderiv ; ++j) {
            ljcs[j] = counter;

            rownnz = 0;

            for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                k = irsA[i];
                col_pos = k / nrow;
                row_pos = k%nrow;
                if (pos[col_pos] != (-1)) {
                    isinflag=0;
                    for (l=0; l < rownnz; ++l) {
                        if (lirs[counter+l]==row_pos){
                            srA[counter+l] += prA[i]*prV[pos[col_pos]];
                            isinflag = 1;
                            break;
                        }
                        if (lirs[counter+l] > row_pos) {
                            for (m=rownnz; m>l; m--) {
                                isinflag=1;
                                srA[counter+m] = srA[counter+m-1];
                                lirs[counter+m] = lirs[counter+m-1];
                            }
                            srA[counter+l] = prA[i]*prV[pos[col_pos]];
                            lirs[counter+l] = row_pos;
                            ++rownnz;
                            isinflag=1;
                            break;
                        }
                    }
                    if (isinflag == 0) {
                        srA[counter+rownnz] = prA[i]*prV[pos[col_pos]];
                        lirs[counter+rownnz] = row_pos;
                        ++rownnz;
                    }
                }
            }
            counter+=rownnz;
        }
        ljcs[nderiv]=counter;

        if (counter > 0) {
            lirs = mxRealloc(lirs, counter * sizeof(*lirs));
            srA = mxRealloc(srA, counter * sizeof(*srA));
        }
        mxFree(pos);
    }

    /* Set Output */
    plhs[0] = mxCreateSparse(nrow,nderiv,counter,mxREAL);
    if (counter>0) {
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
    mxFree(pointer);
}
