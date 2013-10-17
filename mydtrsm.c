#include "mex.h"
void copy(double *Aout, double *A, int M, int N);
double myabs(double n);
void printout(double *A, int M, int N);
void exchangerow(double *A, int N, int i, int j);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mexPrintf("C PROGRAM START\n");
    double *Ain, *A;
    int *pvtin, *pvt;
    int M,N, i, j, k, tmp;
    int maxind;
    double max;
    /*mxget*/
    Ain = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]);
    pvtin = mxGetPr(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,N,mxREAL);
    A = mxGetPr(plhs[0]);
    pvt = mxGetPr(plhs[1]);
    /*copy*/
    copy(A, Ain, N, N);
    copy(pvt, pvtin, 1, N);
    for (i=0;i<N;i++) {
        maxind=i;
        max=myabs(A[i*N+i]);
        for (j=i+1;j<N;j++) {
            if (myabs(A[i*N+j])>max) {
                maxind=j; max=myabs(A[i*N+j]);
            }
        }
        if (max==0) {
            mexPrintf("LUfactoration failed: coefficient matrix is singular");
            return 0; 
        }
        else {
            if (maxind != i) {
                tmp = pvt[i];
                pvt[i] = pvt[maxind];
                pvt[maxind]=tmp;
                exchangerow(A, N, i, maxind);
            }
        }
        for (j=i+1;j<N;j++) {
            A[i*N+j] = A[i*N+j]/A[i*N+i];
            for (k=i+1;k<N;k++) {
                A[k*N+j] = A[k*N+j] - A[i*N+j] * A[k*N+i];
            }
        }
    }
}
void exchangerow(double *A, int N, int i, int j){
    int k;
    double tmp;
    for (k=0;k<N;k++) {
        tmp = A[k*N+i];
        A[k*N+i] = A[k*N+j];
        A[k*N+j] = tmp;
    }
}
void printout(double *A, int M, int N){
    int i,j;
    for (i=0;i<M;i++) {
        for (j=0;j<N;j++) {
            mexPrintf("%f ", A[j*M+i]);
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");
}
void copy(double *Aout, double *A, int M, int N) {
    int i;
    for (i=0;i<M*N;i++) {
        Aout[i]=A[i];
    }
}
double myabs(double n) {
    return (n>0?n:-n);
}
