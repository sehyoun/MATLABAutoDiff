/* Inputs: prhs[0] = [m x n]  matrix A
 *         prhs[1] = [(m*l) x nderiv] matrix of dB'/dx
 * Outputs: plhs[0] = [(n*l) x nderiv] matrix corresponding to (A*dB/dx)
 *
 * Note: To get dB/dx to dB'/dx, you can call dertransp(dB/dx,m) prior to
 *                             calling matvalXmatder
 *
 * Method: This multiplication does on the fly insert-sort. However,
 *
 * by SeHyoun Ahn, Aug 2016
 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    mwIndex n, m, l, nderiv,nnz, counter;

    /* Read In A */
    mwIndex *irsA, *jcsA;
    double *prA;
    m = mxGetN(prhs[0]);
    n = mxGetM(prhs[0]);
    irsA = mxGetIr(prhs[0]);
    jcsA = mxGetJc(prhs[0]);
    prA = mxGetPr(prhs[0]);

    /* Read in dB'/dx */
    mwIndex *irsB, *jcsB;
    double *prB;
    nderiv = mxGetN(prhs[1]);
    l = mxGetM(prhs[1]);
    irsB = mxGetIr(prhs[1]);
    jcsB = mxGetJc(prhs[1]);
    prB = mxGetPr(prhs[1]);
    l = l/m;

    /* Prepare Output Matrix */
    mwIndex *lirs, *ljcs;
    double *lpr;
    nnz = jcsB[nderiv] + 1;

    plhs[0] = mxCreateSparse(n*l, nderiv, nnz, mxREAL);

    lirs = mxGetIr(plhs[0]);
    ljcs = mxGetJc(plhs[0]);
    lpr = mxGetPr(plhs[0]);

    mwIndex i,j, cloc, cir, k,o,p;
    mwSignedIndex pos;
    ljcs[0] = 0;
    counter = 0;

    /* Multiplication Loop */
    for (i=0; i<nderiv; ++i) {
        cloc = counter;
        for (j=jcsB[i]; j<jcsB[i+1]; ++j) {
            pos = irsB[j]/l;
            for (k=jcsA[pos]; k<jcsA[pos+1]; ++k) {
                cir = (irsB[j]%l)*n + irsA[k];
                if (counter == cloc) {
                    if (counter > nnz-1) {
                        nnz = ((nnz+1)*(nderiv+1)/(i+1));
                        lirs = mxRealloc(lirs, nnz*sizeof(*lirs));
                        lpr = mxRealloc(lpr, nnz*sizeof(*lpr));
                    }
                    lirs[counter] = cir;
                    lpr[counter] = prB[j]*prA[k];
                    ++counter;
                }
                else {
                    for (o = counter; o>cloc; --o) {
                        if (lirs[o-1] == cir) {
                            lpr[o-1] += prB[j]*prA[k];
                            break;
                        }
                        else if (lirs[o-1] < cir) {
                            if (counter > nnz-1) {
                                nnz = ((nnz+1)*(nderiv+1)/(i+1));
                                lirs = mxRealloc( lirs, nnz*sizeof(*lirs));
                                lpr = mxRealloc( lpr, nnz*sizeof(*lpr));
                            }
                            for (p=counter; p>o; --p) {
                                lirs[p] = lirs[p-1];
                                lpr[p] = lpr[p-1];
                            }
                            lirs[o] = cir;
                            lpr[o] = prB[j]*prA[k];
                            ++counter;
                            break;
                        }
                        else if (o == cloc+1) {
                            if (counter > nnz-1) {
                                nnz = ((nnz+1)*(nderiv+1)/(i+1));
                                lirs = mxRealloc( lirs, nnz*sizeof(*lirs));
                                lpr = mxRealloc( lpr, nnz*sizeof(*lpr));
                            }
                            for (p=counter; p>cloc; --p) {
                                lirs[p] = lirs[p-1];
                                lpr[p] = lpr[p-1];
                            }
                            lirs[cloc] = cir;
                            lpr[cloc] = prB[j]*prA[k];
                            ++counter;
                        }
                    }
                }
            }
        }
        ljcs[i+1]=counter;
    }


    /* Prepare Output */
    if (counter > 0) {
        lirs = mxRealloc(lirs, counter*sizeof(*lirs));
        lpr = mxRealloc(lpr, counter*sizeof(*lpr));

        mxSetIr(plhs[0], lirs);
        mxSetPr(plhs[0], lpr);
        mxSetNzmax(plhs[0], counter+1);
    }
}
