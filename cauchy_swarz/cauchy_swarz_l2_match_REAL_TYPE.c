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

void cauchy_swarz_l2_match(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L);



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

    plhs[0]=() mxCreateDoubleMatrix(c_r,c_c,mxREAL);

    
    
    cauchy_swarz_l2_match(I,img_r,img_c,box_arr,box_num, w_arr, threshold, reconst_pat_copy,n_p,m_p, (REAL_TYPE*)mxGetPr(plhs[0]));
    /*    mxFree(box_arr_copy);*/
    mxFree(reconst_pat_copy);
    /*    box_corr_given_indices(I,img_r,img_c,box_arr,box_num,w_arr,n_p,m_p,tpi,tpj,ind_num, mxGetPr(plhs[0]));*/
}

void cauchy_swarz_l2_match(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L){

    printf("0\n");
    REAL_TYPE* I = mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
    REAL_TYPE* img_sqrs =  mxMalloc(img_r*img_c*sizeof(REAL_TYPE));
    
    /* construct integral image + zeros pad (for boundary problems) */
    /* zero(I,(img_r+1)*(img_c+1));*/
    integral_image(img,img_r,img_c,I); 



    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;


    REAL_TYPE* wins_sum_sqrs=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* wins_norms=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* C=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);

    
    
    int k;
    for(k=0;k<img_r*img_c;k++){
        img_sqrs[k]=img[k]*img[k];
    }

    int window[4];
    window[0]=1;
    window[1]=n_p;
    window[2]=1;
    window[3]=m_p;
    REAL_TYPE window_w[1];
    window_w[0]=1.0;
        

    
    /*box_corr(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num, const REAL_TYPE* w_arr,int n_p,int m_p, REAL_TYPE* C)*/

    /**********************************************************************/
    /*IMPORTANT remember to zero the output C! , also zero integral image */
    /**********************************************************************/
    box_corr_matlab_indices(img_sqrs,img_r,img_c,window,1,window_w,n_p,m_p, wins_sum_sqrs);
      
    
    matrix_sqrt(wins_sum_sqrs,c_r*c_c, wins_norms);
    REAL_TYPE pat_sum_sqrs =sum_sqrs(reconst_pat,n_p,m_p);

     
   /*void matrix_add(const REAL_TYPE* m1, const REAL_TYPE* m2,REAL_TYPE addScalar, REAL_TYPE mulScalar1, REAL_TYPE mulScalar2, int size, REAL_TYPE* out){*/
    /*
    ## equivalent to:
    ## L= r2 +  wins_sum_sqrs-2*norm(reconst_filter,'fro')*wins_norm; //L is a matrix with the bound for wach window
    */
    /* L is a lower bound to the l2 distance between the pattern and the window */
    /*REAL_TYPE* L=mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);*/

    REAL_TYPE mul2=(-1)*2*sqrt(pat_sum_sqrs);
    matrix_add(wins_sum_sqrs, wins_norms, pat_sum_sqrs , 1.0,mul2, c_r*c_c, L);

    threshold = threshold*n_p*m_p;
        
    int t=0;

    int* i[2];
    int* j[2];
    
     i[0]=mxMalloc(c_r*c_c*sizeof(int));
     j[0]=mxMalloc(c_r*c_c*sizeof(int));

     i[1]=mxMalloc(c_r*c_c*sizeof(int));
     j[1]=mxMalloc(c_r*c_c*sizeof(int));

    int indicator=0;
    
    int indices_num = find(L, c_r,c_c, threshold,i[0],j[0]);
    
    while(t<box_num){
        
                    
        int box[4];
        box[0]=box_arr[t];
        box[1]=box_arr[t+box_num];
        box[2]=box_arr[t+2*box_num];
        box[3]=box_arr[t+3*box_num];
        REAL_TYPE w[1];
        w[0]=w_arr[t];
         
         /* C=C+box_corr_given_indices_c(int_img,box_arr(t,:),w_arr(t),n_p,m_p,int32(i),int32(j));*/
        /*zero(Ctmp,c_r*c_c);*/
        box_corr_given_indices_ij_by_c_ind_box_by_matlab_ind(I,img_r,img_c,box,1,w, n_p, m_p,i[indicator],j[indicator] ,indices_num, C,0);

    
        /* f=zeros(n_p,m_p);
         * f(box_arr(t,1):box_arr(t,2),box_arr(t,3):box_arr(t,4))=w_arr(t);
         * reconst_filter=reconst_filter-f;
         */
        remove_box_matlab_indices(reconst_pat,n_p,m_p,box[0],box[1],box[2],box[3],w_arr[t]);
        REAL_TYPE mul=-2*sqrt(sum_sqrs(reconst_pat,n_p,m_p));

        add_3_matrices_in_indices(wins_sum_sqrs,1.0,wins_norms,mul,C,-2,pat_sum_sqrs,c_r,i[indicator],j[indicator],indices_num,L);
        t++;

        /*  ind=L<threshold;
         *  [i,j]=find(ind);
         */
        indices_num=find_from_indices(L,c_r,c_c,threshold,indices_num,i[indicator],j[indicator],i[1-indicator],j[1-indicator]);
        indicator=1-indicator;
    }  
    
    mxFree(I);
    mxFree(img_sqrs);
    mxFree(wins_sum_sqrs);
    mxFree(wins_norms);
    mxFree(C);
    mxFree(i[0]);
    mxFree(j[0]);
    mxFree(i[1]);
    mxFree(j[1]);
    
}
