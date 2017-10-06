#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>

//declare variables
double *a, *b;
void sumNumbers(double a, double b, double *op);
void computeTraceTerm(double *inp, double *cov, double *op, int dimx, int dimy, int dimz);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Associate inputs
    double *ftMatrix = mxGetPr(prhs[0]); // The feature vectors
    
    // Figure out dimensions
    const mwSize *dims = mxGetDimensions(prhs[0]);
    int numdims = mxGetNumberOfDimensions(prhs[0]);
    int dimy = (int)dims[0];
    int dimx = (int)dims[1];
    int dimz =  (int)dims[2];
    
    double *covMat = mxGetPr(prhs[1]); // The feature vectors
    
    // output
    plhs[0] = mxCreateDoubleMatrix(dimz, 1, mxREAL);
    double *outArray = mxGetPr(plhs[0]);
    
    /* call the computational routine */
    computeTraceTerm(ftMatrix, covMat, outArray, dimx, dimy, dimz);
}

void computeTraceTerm(double *inp, double *cov, double *op, int dimx, int dimy, int dimz){
    double traceVal;
    int nD = dimx;
    
    for(int i=0; i<dimz; i++){
        traceVal = 0;        
        for(int j=0; j<dimy; j++){
            for(int k=0; k<dimx; k++){
                traceVal += (inp[j+nD*(k+nD*i)] * cov[k+nD*j]);
            }
        }
        op[i] = traceVal;
    }
    
    /*std::cout << dimx << std::endl;
    std::cout << dimy << std::endl;
    std::cout << dimz << std::endl;*/
}