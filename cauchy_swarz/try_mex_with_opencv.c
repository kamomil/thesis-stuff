#include "mex.h"
#include "matrix.h"
/*#include "blas.h"*/
#include "box_corr_lib.h"
#include <stdio.h>
#include <math.h>

#define MAX_IMG_SIZE 1000000
#define MAX_BOX_NUM 3000



/* REMEMBER a RxC matrix is represented such as [1, ... R, R+1, ....  2R ....] */
/* That is, C consecutive blocks each of size R, the ith block represent the ith column */

/*void box_corr_given_indices(const REAL_TYPE* I, int img_r,int img_col, const int* box_arr,int box_num, const REAL_TYPE*
  w_arr,int n_p,int m_p,const int* tpi,const int* tpj , int ind_num, REAL_TYPE* C);*/

void cauchy_swarz_before_loop_opt(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L, int opt);
void cauchy_swarz_l2_match(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L, int opt);
void cauchy_swarz_inside_loop_opt(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L,int old);
void cauchy_swarz_inside_loop_opt_x2(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L,int old);


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


    const REAL_TYPE* reconst_pat = (const REAL_TYPE*)mxGetPr(prhs[1]); 
    int n_p = mxGetM(prhs[4]); 
    int m_p = mxGetN(prhs[4]); 
    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    int I_r=img_r+1;
    int I_c=img_c+1;

    plhs[0]= mxCreateNumericMatrix(c_r,c_c,MX_REAL_TYPE ,mxREAL);

    REAL_TYPE* result_arr =  (REAL_TYPE*)mxGetPr(plhs[0]);
  
  //  Mat m(rows, cols, CV_16U, myShortArray);
  Mat img(5, 6, CV_32F,I);
  Mat templ(2, 2, CV_32F,reconst_pat);

  int result_cols = img.cols - templ.cols + 1;
  int result_rows = img.rows - templ.rows + 1;

  Mat result( result_cols, result_rows, CV_32F );

  matchTemplate(img, templ ,  result, CV_TM_SQDIFF);

  //  float* matData = (float*)myMat.data;
  result_arr = (float*)result.data;                   
}
    
/************************************
 ************************************/
void cauchy_swarz_inside_loop_opt(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L,int orig){

    REAL_TYPE* I = mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
    REAL_TYPE* I_sqrs = mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));

    int_img_and_int_img_sqr(img, img_r, img_c ,I,  I_sqrs);
        
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    REAL_TYPE* wins_sum_sqrs=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* wins_norms=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* C=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);

    REAL_TYPE pat_sum_sqrs =sum_sqrs(reconst_pat,n_p,m_p);
    REAL_TYPE mul2=(-1)*2*sqrt(pat_sum_sqrs);

    
    wins_sums(I_sqrs,img_r,img_c, n_p, m_p,wins_sum_sqrs);
    mxFree(I_sqrs);
    matrix_sqrt2(wins_sum_sqrs,c_r*c_c, wins_norms);
    matrix_add2(wins_sum_sqrs, wins_norms, pat_sum_sqrs , 1.0,mul2, c_r*c_c, L);
    
    threshold = threshold*n_p*m_p;
        
    int t=0;

    int* i[2];
    int* j[2];
    
     i[0]=mxMalloc(c_r*c_c*sizeof(int));
     j[0]=mxMalloc(c_r*c_c*sizeof(int));
     i[1]=mxMalloc(c_r*c_c*sizeof(int));
     j[1]=mxMalloc(c_r*c_c*sizeof(int));

    int indicator=0;
    
    int indices_num = find2(L, c_r,c_c, threshold,i[indicator],j[indicator]);
    /*    printf("inside loop opt t=%d, ind=%d thresh=%f L sum sqrs = %f\n",t,indices_num,threshold,sum_sqrs(L,c_r,c_c));*/
    
    zero2(C,c_c*c_r);            
    while(t<box_num){
        
        indices_num = inside_loop3(I,img_r+1,C,L,wins_sum_sqrs,wins_norms,c_r,c_c,box_arr[t]-1,box_arr[t+box_num]-1,box_arr[t+2*box_num]-1,box_arr[t+3*box_num]-1,w_arr[t],pat_sum_sqrs,reconst_pat,n_p, m_p,i[indicator], j[indicator],i[1-indicator],j[1-indicator], indices_num,threshold);
     
        t++;
        indicator=1-indicator;
    }
    
    mxFree(I);
    mxFree(wins_sum_sqrs);
    mxFree(wins_norms);
    mxFree(C);
    mxFree(i[0]);
    mxFree(j[0]);
    mxFree(i[1]);
    mxFree(j[1]);
    
}
    

void cauchy_swarz_inside_loop_opt_x2(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L,int orig){
    
    REAL_TYPE* I = mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
    REAL_TYPE* I_sqrs = mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));

    int_img_and_int_img_sqr(img, img_r, img_c ,I,  I_sqrs);
        
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    REAL_TYPE* wins_sum_sqrs=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* wins_norms=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* C=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);

    REAL_TYPE pat_sum_sqrs =sum_sqrs(reconst_pat,n_p,m_p);
    REAL_TYPE mul2=(-1)*2*sqrt(pat_sum_sqrs);

    L_before_loop(I_sqrs,img_r,img_c,n_p,m_p, mul2,pat_sum_sqrs, wins_sum_sqrs, wins_norms, L);
    mxFree(I_sqrs);
    
    threshold = threshold*n_p*m_p;
        
    int t=0;

    int* i[2];
    int* j[2];
    
     i[0]=mxMalloc(c_r*c_c*sizeof(int));
     j[0]=mxMalloc(c_r*c_c*sizeof(int));
     i[1]=mxMalloc(c_r*c_c*sizeof(int));
     j[1]=mxMalloc(c_r*c_c*sizeof(int));

    int indicator=0;
    int indices_num = find2(L, c_r,c_c, threshold,i[indicator],j[indicator]);
    
    zero2(C,c_c*c_r);
    int I_r=img_r+1;
            
    while(t<box_num){

        indices_num = inside_loop3(I,img_r+1,C,L,wins_sum_sqrs,wins_norms,c_r,c_c,box_arr[t]-1,box_arr[t+box_num]-1,box_arr[t+2*box_num]-1,box_arr[t+3*box_num]-1,w_arr[t],pat_sum_sqrs,reconst_pat,n_p, m_p,i[indicator], j[indicator],i[1-indicator],j[1-indicator], indices_num,threshold);        
        t++;
        indicator=1-indicator;
    }
    
    mxFree(I);
    mxFree(wins_sum_sqrs);
    mxFree(wins_norms);
    mxFree(C);
    mxFree(i[0]);
    mxFree(j[0]);
    mxFree(i[1]);
    mxFree(j[1]);
    
}
