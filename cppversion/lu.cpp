#include<iostream>
#include<stdlib.h>
#include"f2c.h"
#include"clapack.h"
#include<time.h>
#include<cmath>
#define PI 3.1415926
using namespace std;
double *myrand(int m, int n);
void mydgetrf(long int n, double *A, long int *pvt, int optim);
void exchangerow(double *A, long int n, int i, int j);
double mysumdot(double *a, double *b);
double *subv(double *vect, int flag, int start, int end, long int n);
void mydtrsm( char UPLO, long int n, double *LU, double *B, long int *pvt);
double mynorm(double *b1, double *b2, long int n);
double *mydeepcp(double *a, long int m, long int n);
int main(int argc, char** argv)
{
	long int n, info;
	char SIDE = 'L';
	char UPLO = 'L';
	char TRANSA = 'N';
	char DIAG = 'N';
	double ALPHA = 1.0;
	srand(time(NULL));
	clock_t t1, t2;
	if (argc!=3) {
	    cout<<"USAGE: ./lu [DIMENTION] [LEVEL]"<<endl;
	    cout<<"explaination:\t[DIMENTION] : 1000, 2000, etc"<<endl;
	    cout<<"\t\t[LEVEL]:\n\t\t\t0 ==> calculate my non-optimized version and lapack version"<<endl;
	    cout<<"\t\t\t1 ==> calculate my optimized version and lapack version"<<endl;
	    return 1;
	}
	n=atoi(argv[1]);
	int level = atoi(argv[2]);
	if (n<0 || level>1 || level<0) {
	    cout<<"bad arguments"<<endl;
	    return 1;
	}
	cout<<"dimention of A: "<<n<<"x"<<n<<endl;
	double *A, *B, *Abk, *Bbk;
	//get random matrix and copies
	A = myrand(n, n);
	B = myrand(1, n);
	Abk = mydeepcp(A, n, n);
	Bbk = mydeepcp(B, 1, n);
	long int pvt[n];
	//my algorithm, result in Abk
	if (level==0) cout<<"My Algorithm: non-optimized"<<endl;
	else cout<<"My Algorithm: optimized"<<endl;
	t1 = clock();
	/*
	 * I do some optimization with mydgetrf.
	 * The last parameter is the flag: 
	 * 0 for original algorithm
	 * 1 for optimed algorithm
	 */
	mydgetrf(n, A, pvt, level);
	mydtrsm('L', n, A, B, pvt);
	mydtrsm('U', n, A, B, pvt);
	t2 = clock();
	cout<<"\t"<<double(t2-t1)/CLOCKS_PER_SEC<<" sec"<<endl;
	//use lapack, result in Bbk
	cout<<"Lapack Algorithm"<<endl;
	t1=clock();
	long int mm=1;
	dgetrf_(&n, &n, Abk, &n, pvt, &info);
	dlaswp_(&mm, Bbk, &n, &mm, &n, pvt, &mm);
	UPLO='L';DIAG='U';
	dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &n, &mm, &ALPHA, Abk, &n, Bbk, &n);
	UPLO='U';DIAG='N';
	dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &n, &mm, &ALPHA, Abk, &n, Bbk, &n);
	t2=clock();
	cout<<"\t"<<double(t2-t1)/CLOCKS_PER_SEC<<" sec"<<endl;
	//print norm
	double norm = mynorm(B, Bbk, n);
	cout<<"The norm between My alogorithm and the lapack one is: "<<norm<<endl;
	free(A);free(B);free(Abk);free(Bbk);
	return 0;
}
void mydgetrf(long int n, double *A, long int *pvt, int level) {
/* 
* LU factorization using the algorithm from matlab
* done some optimization here
*/
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
		if (level==1) {
			// the optimized code.
			/* calculate the column from A(i,i) to A(n,i) */
			for (int j=i+1;j<n;j++) {
				A[i*n+j] = A[i*n+j]/A[i*n+i];
			}
			/* calculate the row from A(i,i+1) to A(i, n) */
			for (int j=i+1;j<n;j++) {
				for (int k=0;k<i;k++) {
					A[j*n+i] -= A[k*n+i]*A[j*n+k];
				}
			}
			/* calculate the column from A(i,i+1) to A(n, i+1) */
			for (int j=i+1;j<n;j++) {
				for (int k=0;k<i+1;k++) {
					A[(i+1)*n+j] -= A[k*n+j]*A[(i+1)*n+k];
				}
			}
		} else {
			// the original one
			for (int j=i+1;j<n;j++) {
				A[i*n+j] = A[i*n+j]/A[i*n+i];
				for (int k=i+1;k<n;k++) {
					A[k*n+j] = A[k*n+j] - A[i*n+j]*A[k*n+i];
				}
			}
		}
	}
}
void mydtrsm( char UPLO, long int n, double *A, double *B, long int *pvt) {
/* mydtrsm, use 'U' and 'L' to specify */
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
/* return sub vector */
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
/* do both dot and sum */
	double tmp=0;
	for (int i=1;i<a[0]+1;i++) {
		tmp += a[i]*b[i];
	}
	return tmp;
}

void exchangerow(double *A, long int n, int i, int j) {
	/* exchange row i and row j*/
	for (int k=0;k<n;k++) {
		double tmp = A[k*n+i];
		A[k*n+i] = A[k*n+j];
		A[k*n+j] = tmp;
	}
}
double mynorm(double *b1, double *b2, long int n) {
	/* cal norm */
	double tmp=0;
	for (int i=0;i<n;i++) {
		tmp += abs(b1[i]-b2[i])*abs(b1[i]-b2[i]);
	}
	tmp = sqrt(tmp);
	return tmp;
}
double *myrand(int m, int n) {
	/* get random numbers which are normal distribution where mu=0 */
	int i;
	double *tmp = (double*)malloc(m*n*sizeof(double));
	for (i=0;i<m*n;i++) {
		tmp[i] = sqrt(-2*log(rand()/(RAND_MAX+1.0)))*cos(2*PI*rand()/(RAND_MAX+1.0));
	}
	return tmp;
}
double *mydeepcp(double *a, long int m, long int n) {
	/* deep copy */
	double *tmp = (double*)malloc(m*n*sizeof(double));
	for (int i=0;i<m;i++) {
		for (int j=0;j<n;j++) {
			tmp[j*m+i] = a[j*m+i];
		}
	}
	return tmp;
}
