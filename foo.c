#include "mex.h"
#include<math.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mexPrintf("C PROGRAM START\n");
    int a=-1;
    mexPrintf("%d\n", cabs(a));
}
