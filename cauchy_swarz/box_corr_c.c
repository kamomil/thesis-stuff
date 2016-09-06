#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include "box_corr_lib.h"

#define MAX_IMG_SIZE 1000000
#define MAX_BOX_NUM 3000

/* REMEMBER a RxC matrix is represented such as [1, ... R, R+1, ....  2R ....] */
/* That is, C consecutive blocks each of size R, the ith block represent the ith column */



/* prhs - the array of inputs (each cell contain the pointer and size), will */
/* also have pointer to the output */
void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* Pointer to an mxArray of type double */
    const double* img = mxGetPr(prhs[0]); 
    int img_r = mxGetM(prhs[0]); 
    int img_c = mxGetN(prhs[0]); 


    mxClassID category = mxGetClassID(prhs[1]);

    if(category == mxDOUBLE_CLASS)
        mexErrMsgTxt("box_arr should be int (recieved double)");
    
    const int* box_arr = (int *) mxGetData(prhs[1]);
    int box_num = mxGetM(prhs[1]);

    int from_r = box_arr[0]-1;
    int to_r = box_arr[1]-1;
    int from_c = box_arr[2]-1;
    int to_c = box_arr[3]-1;
  
    const double* w_arr =  mxGetPr(prhs[2]);
    int n_p = (int)mxGetScalar(prhs[3]);
    int m_p = (int)mxGetScalar(prhs[4]); 

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    
    plhs[0]=mxCreateDoubleMatrix(c_r,c_c,mxREAL);
    box_corr_matlab_indices(img,img_r,img_c, box_arr,box_num, w_arr,n_p,m_p,  mxGetPr(plhs[0]));

}


    
/*
void print_matrix(double* m,int r,int c){

    int i;
    int j;
    int k;
    for(i=0;i<r;i++){
        j=0;
        for(k=0;k<c;k++){
            printf("%f ",m[i+j]);
            j+=r;
        }
        printf("\n");
    }
}
*/

