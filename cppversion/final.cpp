#include<math.h>
#include<iostream>
#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include<time.h>
#include<cmath>
#define PI 3.1415926
#define rd (rand()/(RAND_MAX+1.0))
using namespace std;
double *myrand(int m, int n);
void printout(const char* str, double *v, int m, int n);
void iprintout(const char* str, long int *v, int m, int n);
void lihebi(double *A, double *B, double *x, long int n);
double *getL(double* A, long int n);
double *getU(double* A, long int n);
double singlerand();
void getLU(double *A, double *L, double *U, long int n);
int mydgetrf(long int n, double *A, long int *pvt);
void exchangerow(double *A, long int n, int i, int j);
double mysumdot(double *a, double *b);
double *subv(double *vect, int flag, int start, int end, long int n);
void mydtrsm( char UPLO, long int n, double *LU, double *B, long int *pvt);
double mynorm(double *b1, double *b2, long int n);
int main()
{
	double tmp;
	long int n, i, j, k, info;
	//double *A, *B;
	double dtmp;
	char SIDE = 'L';
	char UPLO = 'L';
	char TRANSA = 'N';
	char DIAG = 'N';
	double ALPHA = 1.0;
	char flag='N';
	srand(time(NULL));
	n = 4;
	double *A = (double*)malloc(n*n*sizeof(double));
	double *B = (double*)malloc(n*sizeof(double));
	for (int i=0;i<n*n;i++) {
	    A[i] = singlerand();
	}
	for (int i=0;i<n;i++) {
	    B[i] = singlerand();
	}
	//A[0]=2;A[1]=1;A[2]=3;A[3]=2;B[0]=2;B[1]=3;
	//A[0]=1;A[1]=2;A[2]=2;A[3]=3;B[0]=3;B[1]=2;
	//for (int i=0;i<9;i++) A[i]=i+1; for (int i=0;i<3;i++) B[i]=i+1;
	double *Abk = (double*)malloc(n*n*sizeof(double));
	double *Bbk = (double*)malloc(n*sizeof(double));
	for (int i=0;i<n;i++) {
	    Bbk[i] = B[i];
	    for (int j=0;j<n;j++) {
		Abk[j*n+i] = A[j*n+i];
	    }
	}
	long int pvt[n];
	printout("originA=======", A, n, n);
	printout("originB=======", B, 1, n);
	//printout("Abk =========", Abk, n, n);
//#define MY
#ifdef MY
	cout<<"my"<<endl;
	mydgetrf(n, A, pvt);
	mydtrsm('L', n, A, B, pvt);
	mydtrsm('U', n, A, B, pvt);
	//printout("B===========", B, 1, n);
#else
	cout<<"lapack"<<endl;
	long int mm=1;
//#define GESV
#ifdef GESV
	dgesv_(&n, &mm, A, &n, pvt, B, &n, &info);
#else
	dgetrf_(&n, &n, A, &n, pvt, &info);
	dlaswp_(&mm, B, &n, &mm, &n, pvt, &mm);
	UPLO='L';
	DIAG='U';
	dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &n, &mm, &ALPHA, A, &n, B, &n);
	UPLO='U';
	DIAG='N';
	dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &n, &mm, &ALPHA, A, &n, B, &n);
#endif
#endif
	//printout("Abk =========", Abk, n, n);
	iprintout("pvt =======", pvt, 1, n);
	lihebi(Abk, Bbk, B, n);
	//double norm = mynorm(Bbk, Bbk, n);
	//cout<<norm<<endl;
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
int mydgetrf(long int n, double *A, long int *pvt) {
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
	    return 0;
	}
	else {
	    if (maxind != 1) {
		long int tmp = pvt[i];
		pvt[i] = pvt[maxind];
		pvt[maxind] = tmp;
		exchangerow(A, n, i, maxind);
	    }
	}
#define optim
#ifdef optim
	for (int j=i+1;j<n;j++) {
	    for (int k=0;k<i;k++) {
	    	A[i*n+j] -= A[k*n+j]*A[i*n+k];
	    }
	    A[i*n+j] = A[i*n+j]/A[i*n+i];
	}
	for (int j=i+1;j<n;j++) {
	    for (int k=0;k<i+1;k++) {
		A[j*n+i+1] -= A[k*n+i+1]*A[j*n+k];
	    }
	}
#else
	for (int j=i+1;j<n;j++) {
	    A[i*n+j] = A[i*n+j]/A[i*n+i];
	    for (int k=i+1;k<n;k++) {
		A[k*n+j] = A[k*n+j] - A[i*n+j]*A[k*n+i];
	    }
	}
#endif
    }
}
void exchangerow(double *A, long int n, int i, int j) {
    for (int k=0;k<n;k++) {
	double tmp = A[k*n+i];
	A[k*n+i] = A[k*n+j];
	A[k*n+j] = tmp;
    }
}
void lihebi(double *A, double *B, double *x, long int n) {
    double tmp;
    for (int i=0;i<n;i++) {
	tmp = 0;
	for (int j=0;j<n;j++) {
	    cout<<"\t"<<A[j*n+i]<<" "<<x[j];
	    tmp += A[j*n+i]*x[j];
	}
	cout<<"==>>> "<<tmp<<endl;
    }
}
double mynorm(double *b1, double *b2, long int n) {
    double tmp=0;
    for (int i=0;i<n;i++) {
	tmp += abs(b1[i]-b2[i]);
    }
    return tmp;
}
void getLU(double *A, double *L, double *U, long int n) {
    for (int i=0;i<n;i++) {
	for (int j=0;j<n;j++) {
	    if (i<j) {
		L[j*n+i] = 0;
		U[j*n+i] = A[j*n+i];
	    }else if (i==j) {
		L[j*n+i] = 1;
		U[j*n+i] = A[j*n+i];
	    } else {
		L[j*n+i] = A[j*n+i];
		U[j*n+i] = 0;
	    }
	}
    }
}
double *getL(double* A, long int n) {
    double *L = (double*)malloc(n*n*sizeof(double));
    for (int i=0;i<n;i++) {
	for (int j=0;j<n;j++) {
	    if (i<j) {
		L[j*n+i] = 0;
	    }else if (i==j) {
		L[j*n+i] = 1;
	    } else {
		L[j*n+i] = A[j*n+i];
	    }
	}
    }
    return L;
}
double *getU(double* A, long int n) {
    double *L = (double*)malloc(n*n*sizeof(double));
    for (int i=0;i<n;i++) {
	for (int j=0;j<n;j++) {
	    if (i>j) {
		L[j*n+i] = 0;
	    } else {
		L[j*n+i] = A[j*n+i];
	    }
	}
    }
    return L;
}
void printout(const char* str, double *v, int m, int n) {
	int i,j;
	cout<<str<<endl;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			cout<<v[j*m+i]<<"\t";
		}
		cout<<endl;
	}
}
void iprintout(const char* str, long int *v, int m, int n) {
	int i,j;
	cout<<str<<endl;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			cout<<v[j*m+i]<<"\t";
		}
		cout<<endl;
	}
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
