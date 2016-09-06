#include "mex.h"
#include "matrix.h"
/*#include "blas.h"*/
#include "../cauchy_lib/cauchy_swarz_l2_match.h"
#include <stdio.h>
#include <math.h>

/* REMEMBER a RxC matrix is represented such as [1, ... R, R+1, ....  2R ....] */
/* That is, C consecutive blocks each of size R, the ith block represent the ith column */

/* prhs - the array of inputs (each cell contain the pointer and size), will */
/* also have pointer to the output */
void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* Pointer to an mxArray of type REAL_TYPE */
    const REAL_TYPE* I = (const REAL_TYPE*)mxGetPr(prhs[0]); 
    int img_r = mxGetM(prhs[0]); 
    int img_c = mxGetN(prhs[0]); 
   
    if(mxGetClassID(prhs[1]) == mxDOUBLE_CLASS)
        mexErrMsgTxt("box_arr should be int (recieved REAL_TYPE)");
    
    const int* box_arr = (int *) mxGetData(prhs[1]);
    int box_num = mxGetM(prhs[1]);
    
    /*    int* box_arr_copy = mxMalloc(sizeof(int)*4*box_num);
    int k=0;
    for(k=0;k<4*box_num;k++){
        box_arr_copy[k]=box_arr[k];
        }*

        indices_mat2c(box_arr_copy,4*box_num); */
    const REAL_TYPE* w_arr =  (const REAL_TYPE*)mxGetPr(prhs[2]);    
    REAL_TYPE threshold = mxGetScalar(prhs[3]);
    const REAL_TYPE* reconst_pat = (const REAL_TYPE*)mxGetPr(prhs[4]);
    
    int n_p = mxGetM(prhs[4]); 
    int m_p = mxGetN(prhs[4]); 

    REAL_TYPE* reconst_pat_copy = mxMalloc(sizeof(REAL_TYPE)*n_p*m_p);    
    int k;
    for(k=0;k<n_p*m_p;k++){
        reconst_pat_copy[k]=reconst_pat[k];
    }
    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int I_r=img_r+1;
    int I_c=img_c+1;

    plhs[0]= mxCreateNumericMatrix(c_r,c_c,MX_REAL_TYPE ,mxREAL);
    cauchy_swarz_inside_loop_opt(I,img_r,img_c,box_arr,box_num, w_arr, threshold, reconst_pat_copy,n_p,m_p, (REAL_TYPE*)mxGetPr(plhs[0]));
    
    mxFree(reconst_pat_copy);
}
