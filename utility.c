#include "mex.h"
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

mysumdot(double *a, double *b) {
    double tmp=0;
    int i;
    for (i=1;i<a[0]+1;i++) {
        tmp+=a[i]*b[i];
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
