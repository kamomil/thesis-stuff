#include "cauchy_swarz_ncc_match.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define SEC_RES 1000

int cauchy_swarz_ncc_clean(const REAL_TYPE* I,const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p, const REAL_TYPE* norms_of_wins_minus_mu,const REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* i_loc,int* j_loc, REAL_TYPE* vals, REAL_TYPE* U){

   
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    /*
#ifdef MEX
    REAL_TYPE* I =  mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
#else
    REAL_TYPE* I =  malloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));        
#endif
    integral_image_opt(img,img_r, img_c , I);
    */
    
    int indexes_num = U_before_loop_and_find(I,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,i_loc,j_loc,threshold*residual_pat_norms[0]-residual_pat_norms[1]);

    box_arr=box_arr+4;
    int total_indexes_num=c_r*c_c;
    REAL_TYPE cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);

    int t=1;
    while(t<box_num && cur_frac>frac_for_direct){
        int f_r=(*box_arr)-1;
        int t_r=*(box_arr+1)-1;
        int f_c=*(box_arr+2)-1;
        int t_c=*(box_arr+3)-1;
        box_arr=box_arr+4;
        /*prev_indexes_num=indexes_num;*/
        indexes_num = ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,i_loc,j_loc,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1]);
        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        t++;  
    }
    direct_ncc(img,img_r,img_c,pat,n_p,m_p, i_loc,j_loc,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);

    /*
#ifdef MEX
     mxFree(I);
#else
     free(I);
#endif
    */
     return indexes_num;
}
 
