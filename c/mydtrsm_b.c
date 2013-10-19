#include "mex.h"
void printout(double *A, int M, int N);
double mysumdot(double *a, double *b);
double* subv(double *vect, int flag, int start, int end, int N);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mexPrintf("C PROGRAM START\n");
    double *A, *y, *x;
    int i, j, tmp, N;
    double *a,*b;
    /*para in: A, y*/
    /*para out: x*/
    /*mxget*/
    A = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]);
    y = mxGetPr(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
    x = mxGetPr(plhs[0]);
    /*back*/
    x[N-1] = y[N-1] / A[(N-1)*N+N-1];
    for (i=N-2;i>=0;i--) {
        a = subv(x, -1, i+1, N-1, N);
        b = subv(A, i, i+1, N-1, N);
        x[i] = (y[i] - mysumdot(a, b)) /A[i*N+i];
        free(a);
        free(b);
    }
}
double newsubv(double *vect, int flag, int start, int end, int N) {
    int i;
    double tmp;
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
    /*flag is row number*/
    int i;
    /*double tmp[end-start+2];*/
    double *tmp = malloc((end-start+2)*sizeof(double));
    tmp[0]=end-start+1;
    if (flag==-1) {
        /*1 dimension*/
        for (i=0;i<end-start+1;i++) {
            tmp[i+1] = vect[start+i];
        }
    }
    else {
        /*matrix*/
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
