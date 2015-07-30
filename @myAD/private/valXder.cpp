/**************************************************************\
*	Columnwise array multiplication between Matrix and Vector  *
\**************************************************************/

#include "mex.h"
#define	Out	plhs[0]
#define	Vec	prhs[0]
#define	Mat	prhs[1]

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{	
	if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
	switch(nrhs) {
    default: {
    	mexErrMsgTxt("Incorrect number of arguments.");
        }
    case 0: {
        mexPrintf("\nColumnwise dot-product matrix with vector.\n\n");
        break;
        }
	case 2: {
		int m=mxGetN(Vec);
		int n=mxGetN(Mat);
		Out=mxCreateDoubleMatrix(0,0,mxREAL);
		if (m == 1 && n > 0) {
            int i, j, k=0;
			const int* e=mxGetDimensions(Mat);
			double* or;
            double* inMat=mxGetPr(Mat);
            double* inVec=mxGetPr(Vec);
            mxSetDimensions(Out,e,2);
    		mxSetPr(Out,or=(double*) mxMalloc(e[0]*e[1]*sizeof(double)));
            for (i=0; i<e[1]; i++) {
                for (j=0; j<e[0]; j++) {
                    or[k] = inMat[k]*inVec[j];
                    k++;
                }
            }
		}
        }
	}
}
