#include "mex.h"
#include "matrix.h"
/*#include "blas.h"*/
#include "../cauchy_lib/cauchy_swarz_ncc_match.h"
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
    const REAL_TYPE* img = (const REAL_TYPE*)mxGetPr(prhs[0]); 
    int img_r = mxGetM(prhs[0]); 
    int img_c = mxGetN(prhs[0]); 
   
    if(mxGetClassID(prhs[1]) == mxDOUBLE_CLASS)
        mexErrMsgTxt("box_arr should be int (recieved REAL_TYPE)");
    
    const int* box_arr = (int *) mxGetData(prhs[1]);
    int box_num = mxGetN(prhs[1]);
    

    const REAL_TYPE* w_arr =  (const REAL_TYPE*) mxGetPr(prhs[2]);    
    /*    REAL_TYPE threshold = (REAL_TYPE) mxGetScalar(prhs[3]);*/
    /*double  thresholdd = mxGetScalar(prhs[3]);*/
    const REAL_TYPE* th_ptr = (const REAL_TYPE*)mxGetPr(prhs[3]);
    REAL_TYPE threshold = th_ptr[0];
    
    const REAL_TYPE* pat = (const REAL_TYPE*)mxGetPr(prhs[4]);
    const REAL_TYPE* norms_of_wins_minus_mu = (const REAL_TYPE*)mxGetPr(prhs[5]);
    const REAL_TYPE* wins_mu = (const REAL_TYPE*)mxGetPr(prhs[6]);
    const REAL_TYPE* residual_pat_norms = (const REAL_TYPE*)mxGetPr(prhs[7]);
    /*REAL_TYPE frac_for_direct = (REAL_TYPE )mxGetScalar(prhs[8]);*/
    const REAL_TYPE* f_ptr = (const REAL_TYPE*)mxGetPr(prhs[8]);
    REAL_TYPE frac_for_direct = f_ptr[0];

    /*
    REAL_TYPE threshold = (REAL_TYPE) thresholdd;
    */
    /*    printf("thd %f frac %f\n",threshold,frac_for_direct);*/
    
    int n_p = mxGetM(prhs[4]); 
    int m_p = mxGetN(prhs[4]); 

    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int I_r=img_r+1;
    int I_c=img_c+1;

    plhs[0]= mxCreateNumericMatrix(c_r,c_c,MX_REAL_TYPE ,mxREAL);
    plhs[1]= mxCreateNumericMatrix(c_r*c_c,1,mxINT32_CLASS ,mxREAL);
    plhs[2]= mxCreateNumericMatrix(c_r*c_c,1,mxINT32_CLASS ,mxREAL);
    plhs[3]= mxCreateNumericMatrix(c_r*c_c,1,MX_REAL_TYPE ,mxREAL);
    
    plhs[4]= mxCreateNumericMatrix(1,box_num, mxINT32_CLASS ,mxREAL);
    plhs[5]= mxCreateNumericMatrix(1,box_num, MX_REAL_TYPE ,mxREAL);
    

    /*    
    int ind_num=cauchy_swarz_ncc(const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p, const REAL_TYPE* norms_of_wins_minus_mu,const REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,(int*)mxGetPr(plhs[1]),(int*)mxGetPr(plhs[2]),(REAL_TYPE*)mxGetPr(plhs[3]), (REAL_TYPE*)mxGetPr(plhs[0]));
    */
    
    int ind_num=cauchy_swarz_ncc(img,img_r,img_c,pat,n_p,m_p,norms_of_wins_minus_mu,wins_mu, box_arr, w_arr, box_num, threshold, residual_pat_norms,frac_for_direct,(int*)mxGetPr(plhs[1]),(int*)mxGetPr(plhs[2]),(REAL_TYPE*)mxGetPr(plhs[3]), (REAL_TYPE*)mxGetPr(plhs[0]),(int*)mxGetPr(plhs[4]),(REAL_TYPE*)mxGetPr(plhs[5]));


    plhs[6] = mxCreateDoubleScalar(ind_num);
    /*    printf("returning from mexFunction %d\n",ind_num);*/
}
