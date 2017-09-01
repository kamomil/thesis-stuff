#include "cauchy_swarz_ncc_match.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define SEC_RES 1000



int cauchy_swarz_ncc(const REAL_TYPE* I,const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p, const REAL_TYPE* norms_of_wins_minus_mu,const REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* i_loc,int* j_loc, REAL_TYPE* vals, REAL_TYPE* U, int* ind_frac_iter, double* iter_run_time,double* dconv_run_time){

    struct timeval  first_loop_tv1, first_loop_tv2;
    struct timeval  loop_tv1, loop_tv2;
    struct timeval  dconv_tv1, dconv_tv2;
    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    
    gettimeofday(&first_loop_tv1,NULL);
    int indexes_num = U_before_loop_and_find(I,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,i_loc,j_loc,threshold*residual_pat_norms[0]-residual_pat_norms[1]);
    gettimeofday(&first_loop_tv2,NULL);
    iter_run_time[0]=SEC_RES * ((double) (first_loop_tv2.tv_usec - first_loop_tv1.tv_usec) / 1000000 +(double) (first_loop_tv2.tv_sec - first_loop_tv1.tv_sec));

    /*printf("1, run time=%f\n",iter_run_time[0]);*/
    ind_frac_iter[0]=indexes_num;
    box_arr=box_arr+4;
        int total_indexes_num=c_r*c_c;
    REAL_TYPE cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);

    int t=1;
    /*    printf("t=%d ind num=%d  cur_frac=%f  \n",t,indexes_num,cur_frac);*/
    while(t<box_num && cur_frac>frac_for_direct){
        gettimeofday(&loop_tv1,NULL);
        int f_r=(*box_arr)-1;
        int t_r=*(box_arr+1)-1;
        int f_c=*(box_arr+2)-1;
        int t_c=*(box_arr+3)-1;
        box_arr=box_arr+4;
        /*prev_indexes_num=indexes_num;*/
        indexes_num = ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,i_loc,j_loc,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1]);
        ind_frac_iter[t]=indexes_num;
        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        gettimeofday(&loop_tv2,NULL);

        /*        long elapsed = (loop_tv2.tv_sec-loop_tv1.tv_sec)*SEC_RES + loop_tv2.tv_usec-loop_tv1.tv_usec;*/
        iter_run_time[t]=SEC_RES * ((double) (loop_tv2.tv_usec - loop_tv1.tv_usec) / 1000000 +(double) (loop_tv2.tv_sec - loop_tv1.tv_sec));

        /*
        printf ("Total Box time = %f seconds, \n",(double) (loop_tv2.tv_usec - loop_tv1.tv_usec) / SEC_RES +(double) (loop_tv2.tv_sec - loop_tv1.tv_sec));
        printf ("Total Box time = %f seconds, \n",(REAL_TYPE) (loop_tv2.tv_usec - loop_tv1.tv_usec) / SEC_RES +(REAL_TYPE) (loop_tv2.tv_sec - loop_tv1.tv_sec));
        */

        /*printf("%d, run time=%f\n",t,iter_run_time[t]);*/
        /*printf("t=%d ind num=%d  cur_frac=%f  \n",t,indexes_num,cur_frac);*/
        t++;
        
    }

    gettimeofday(&dconv_tv1,NULL);
    direct_ncc(img,img_r,img_c,pat,n_p,m_p, i_loc,j_loc,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);
    gettimeofday(&dconv_tv2,NULL);
    *dconv_run_time= SEC_RES * ((double) (dconv_tv2.tv_usec - dconv_tv1.tv_usec) / 1000000 +(double) (dconv_tv2.tv_sec - dconv_tv1.tv_sec));

    if(*dconv_run_time <0){
        printf("dconv run time=%f\n",*dconv_run_time);
    }
    /*
    for(t=0;t<box_num;t++){
        printf("%f ",iter_run_time[t]);
    }
    printf("\n");
    */

    /*
#ifdef MEX
     mxFree(I);
#else
     free(I);
#endif
    */
     return indexes_num;
}
