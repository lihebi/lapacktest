#include<math.h>
#include<iostream>
#include<stdlib.h>
#include "f2c.h"
#include"clapack.h"
#include<time.h>
#include<cmath>
#define PI 3.1415926
#define rd (rand()/(RAND_MAX+1.0))
using namespace std;
double *myrand(int m, int n);
double singlerand();
void mydgetrf(long int n, double *A, long int *pvt, int optim);
void exchangerow(double *A, long int n, int i, int j);
double mysumdot(double *a, double *b);
double *subv(double *vect, int flag, int start, int end, long int n);
void mydtrsm( char UPLO, long int n, double *LU, double *B, long int *pvt);
double mynorm(double *b1, double *b2, long int n);
int main()
{
	long int n, info;
	char SIDE = 'L';
	char UPLO = 'L';
	char TRANSA = 'N';
	char DIAG = 'N';
	double ALPHA = 1.0;
	srand(time(NULL));
	n = 2000;
	double *A = (double*)malloc(n*n*sizeof(double));
	double *B = (double*)malloc(n*sizeof(double));
	for (int i=0;i<n*n;i++) {
	    A[i] = singlerand();
	}
	for (int i=0;i<n;i++) {
	    B[i] = singlerand();
	}
	double *Abk = (double*)malloc(n*n*sizeof(double));
	double *Bbk = (double*)malloc(n*sizeof(double));
	for (int i=0;i<n;i++) {
	    Bbk[i] = B[i];
	    for (int j=0;j<n;j++) {
		Abk[j*n+i] = A[j*n+i];
	    }
	}
	long int pvt[n];
//---------------------------------
	cout<<"my"<<endl;
	mydgetrf(n, A, pvt, 1);
	mydtrsm('L', n, A, B, pvt);
	mydtrsm('U', n, A, B, pvt);
//--------------------------------
	cout<<"lapack"<<endl;
	long int mm=1;
	//dgesv_(&n, &mm, A, &n, pvt, B, &n, &info);
	dgetrf_(&n, &n, Abk, &n, pvt, &info);
	dlaswp_(&mm, Bbk, &n, &mm, &n, pvt, &mm);
	UPLO='L';
	DIAG='U';
	dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &n, &mm, &ALPHA, Abk, &n, Bbk, &n);
	UPLO='U';
	DIAG='N';
	dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &n, &mm, &ALPHA, Abk, &n, Bbk, &n);
//---------------------------------
	double norm = mynorm(B, Bbk, n);
	cout<<norm<<endl;
	return 0;
}
void mydtrsm( char UPLO, long int n, double *A, double *B, long int *pvt) {
    double y[n];
    double *va, *vb;
    if (UPLO=='L') {
    	y[0] = B[pvt[0]];
    	for (int i=1;i<n;i++) {
		va = subv(y, -1, 0, i-1, n);
		vb = subv(A, i, 0, i-1, n);
		y[i] = B[pvt[i]] - mysumdot(va, vb);
		free(va);
		free(vb);
    	}
	for (int i=0;i<n;i++) {
	    B[i] = y[i];
	}
    } else if (UPLO=='U') {
	y[n-1] = B[n-1] / A[(n-1)*n+n-1];
	for (int i=n-2;i>=0;i--) {
	    va = subv(y, -1, i+1, n-1, n);
	    vb = subv(A, i, i+1, n-1, n);
	    y[i] = (B[i] - mysumdot(va, vb)) / A[i*n+i];
	    free(va);
	    free(vb);
	}
	for (int i=0;i<n;i++) {
	    B[i] = y[i];
	}
    }
}
double *subv(double *vect, int flag, int start, int end, long int n) {
    double *tmp = (double*)malloc((end-start+2)*sizeof(double));
    tmp[0]=end-start+1;
    if (flag==-1) {
	for (int i=0;i<end-start+1;i++) {
	    tmp[i+1] = vect[start+i];
	}
    } else {
	for (int i=0;i<end-start+1;i++) {
	    tmp[i+1]=vect[(start+i)*n+flag];
	}
    }
    return tmp;
}
double mysumdot(double *a, double *b) {
    double tmp=0;
    for (int i=1;i<a[0]+1;i++) {
	tmp += a[i]*b[i];
    }
    return tmp;
}
void mydgetrf(long int n, double *A, long int *pvt, int optim) {
    int maxind;
    for (int i=0;i<n;i++) {
	pvt[i]=i;
    }
    double max;
    for (int i=0;i<n;i++) {
	maxind = i;
	max = abs(A[i*n+i]);
	for (int j=i+1;j<n;j++) {
	    if (abs(A[i*n+j])>max) {
		maxind=j; max=abs(A[i*n+j]);
	    }
	}
	if (max==0) {
	    cout<<"LUfactoration failed: coefficient matrix is singular"<<endl;
	    return;
	}
	else {
	    if (maxind != 1) {
		long int tmp = pvt[i];
		pvt[i] = pvt[maxind];
		pvt[maxind] = tmp;
		exchangerow(A, n, i, maxind);
	    }
	}
	if (optim==1) {
		for (int j=i+1;j<n;j++) {
		    A[i*n+j] = A[i*n+j]/A[i*n+i];
		}
		for (int j=i+1;j<n;j++) {
		    for (int k=0;k<i;k++) {
			A[j*n+i] -= A[k*n+i]*A[j*n+k];
		    }
		}
		for (int j=i+1;j<n;j++) {
		    for (int k=0;k<i+1;k++) {
			A[(i+1)*n+j] -= A[k*n+j]*A[(i+1)*n+k];
		    }
		}
	} else {
		for (int j=i+1;j<n;j++) {
		    A[i*n+j] = A[i*n+j]/A[i*n+i];
		    for (int k=i+1;k<n;k++) {
			A[k*n+j] = A[k*n+j] - A[i*n+j]*A[k*n+i];
		    }
		}
	}
    }
}
void exchangerow(double *A, long int n, int i, int j) {
    for (int k=0;k<n;k++) {
	double tmp = A[k*n+i];
	A[k*n+i] = A[k*n+j];
	A[k*n+j] = tmp;
    }
}
double mynorm(double *b1, double *b2, long int n) {
    double tmp=0;
    for (int i=0;i<n;i++) {
	tmp += abs(b1[i]-b2[i]);
    }
    return tmp;
}
double *myrand(int m, int n) {
	int i;
	double *tmp = (double*)malloc(m*n*sizeof(double));
	for (i=0;i<m*n;i++) {
		tmp[i] = sqrt(-2*log(rand()/(RAND_MAX+1.0)))*cos(2*PI*rand()/(RAND_MAX+1.0));
	}
	return tmp;
}
double singlerand() {
	return sqrt(-2*log(rand()/(RAND_MAX+1.0)))*cos(2*PI*rand()/(RAND_MAX+1.0));
}
