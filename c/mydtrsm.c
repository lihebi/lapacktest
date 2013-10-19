#include "mex.h"
void printout(double *A, int M, int N);
double* subv(double *vect, int flag, int start, int end, int N);
double mysumdot(double *a, double *b);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mexPrintf("C PROGRAM START\n");
    double *A, *b, *y, *x;
    double *pvt, *va, *vb;
    int i, tmp;
    int flag;
    int N;
    /*para in: A, b, pvt*/
    /*para out: y*/
    /*mxget*/
    A = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]);
    b = mxGetPr(prhs[1]);
    pvt= mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
    y = mxGetPr(mxCreateDoubleMatrix(1,N,mxREAL));
    x = mxGetPr(plhs[0]);
    /*forward*/
    tmp = pvt[0]-1;
    y[0] = b[tmp];
    for (i=1;i<N;i++) {
        tmp=pvt[i]-1;
        va = subv(y, -1, 0, i-1, N);
        vb = subv(A, i, 0, i-1, N);
        y[i] = b[tmp] - mysumdot(va, vb);
        free(va);
        free(vb);
    }
    /*back*/
    x[N-1] = y[N-1] / A[(N-1)*N+N-1];
    for (i=N-2;i>=0;i--) {
        va = subv(x, -1, i+1, N-1, N);
        vb = subv(A, i, i+1, N-1, N);
        x[i] = (y[i] - mysumdot(va, vb)) /A[i*N+i];
        free(va);
        free(vb);
    }
}
double mysumdot(double *a, double *b) {
    double tmp=0;
    int i;
    for (i=1;i<a[0]+1;i++) {
        tmp += a[i]*b[i];
    }
    return tmp;
}
double* subv(double *vect, int flag, int start, int end, int N) {
    int i;
    double *tmp = malloc((end-start+2)*sizeof(double));
    tmp[0]=end-start+1;
    if (flag==-1) {
        for (i=0;i<end-start+1;i++) {
            tmp[i+1] = vect[start+i];
        }
    }
    else {
        for (i=0;i<end-start+1;i++) {
            tmp[i+1]=vect[(start+i)*N+flag];
        }
    }
    return tmp;
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
