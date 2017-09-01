#include "cauchy_swarz_ncc_match.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define SEC_RES 1000

int cauchy_swarz_ncc_linear_idx(const REAL_TYPE* I,const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p, const REAL_TYPE* norms_of_wins_minus_mu,const REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* lidx, REAL_TYPE* vals, REAL_TYPE* U, int* ind_frac_iter, double* iter_run_time,double* dconv_run_time){

    struct timeval  first_loop_tv1, first_loop_tv2;
    struct timeval  loop_tv1, loop_tv2;
    struct timeval  dconv_tv1, dconv_tv2;
    
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
    
    gettimeofday(&first_loop_tv1,NULL);
    int indexes_num = linear_idx_U_before_loop_and_find(I,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,lidx,threshold*residual_pat_norms[0]-residual_pat_norms[1]);
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
   
        indexes_num = lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1]);       
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
    lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);
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

/*********************************************************************************************************/


int cauchy_swarz_ncc(const REAL_TYPE* I,const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p, const REAL_TYPE* norms_of_wins_minus_mu,const REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* i_loc,int* j_loc, REAL_TYPE* vals, REAL_TYPE* U, int* ind_frac_iter, double* iter_run_time,double* dconv_run_time){

    struct timeval  first_loop_tv1, first_loop_tv2;
    struct timeval  loop_tv1, loop_tv2;
    struct timeval  dconv_tv1, dconv_tv2;
    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    /*
#ifdef MEX
    REAL_TYPE* I =  mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
#else
    REAL_TYPE* I =/ malloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));        
#endif
    integral_image_opt(img,img_r, img_c , I);
    */
    
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

/*********************************************************************************************************/

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
    
/*************************************************************************************************************************/



/*box arr is [f_r0,t_r0,f_c0,t_c0,f_r1,t_r1,f_c1,t_c1,f_r2,t_r2,f_c2,t_c2,...]*/
int cauchy_swarz_ncc_linear_idx_clean(const REAL_TYPE* I,const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p, const REAL_TYPE* norms_of_wins_minus_mu,const REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* lidx, REAL_TYPE* vals, REAL_TYPE* U){

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    /*    printf("\n%d %d %d %d %f\n",box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0]);*/
    int indexes_num = linear_idx_U_before_loop_and_find(I,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,lidx,threshold*residual_pat_norms[0]-residual_pat_norms[1]);
    
    box_arr=box_arr+4;
    int total_indexes_num=c_r*c_c;
    REAL_TYPE cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
    /*printf("cur frac= %f thresh=%f\n",cur_frac,threshold*residual_pat_norms[0]-residual_pat_norms[1]);*/
    
    int t=1;

    /*    printf("before while %f %d \n",threshold,box_num);*/
    /*    printf("before while %f %d %f\n",cur_frac,indexes_num,frac_for_direct);*/
    while(t<box_num && cur_frac>frac_for_direct){
        int f_r=(*box_arr)-1;
        int t_r=*(box_arr+1)-1;
        int f_c=*(box_arr+2)-1;
        int t_c=*(box_arr+3)-1;
        box_arr=box_arr+4;
        /*prev_indexes_num=indexes_num;*/
        indexes_num = lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1]);          
        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        t++;
        printf("t=%d cur_frac=%f\n",t,cur_frac);
    }
    
    /*printf("threshold=%f frac_for_direct=%f t=%d cur_frac=%f indexes_num=%d\n",threshold,frac_for_direct,t,cur_frac,indexes_num);*/
    /*    printf("t=%d cur_frac=%f indexes_num=%d\n",t,cur_frac,indexes_num);*/
    lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);

    /*
#ifdef MEX
     mxFree(I);
#else
     free(I);
#endif
    */
     return indexes_num;
}

/*************************************************************************************************************************/



/*box arr is [f_r0,t_r0,f_c0,t_c0,f_r1,t_r1,f_c1,t_c1,f_r2,t_r2,f_c2,t_c2,...]*/
int cauchy_swarz_ncc_linear_idx_clean2(const REAL_TYPE* I,const REAL_TYPE* I2, const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p,REAL_TYPE* norms_of_wins_minus_mu,REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* lidx, REAL_TYPE* vals, REAL_TYPE* U,int* iters_num_ptr){

   printf("XXX cauchy_swarz_ncc_linear_idx_clean2\n");
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;   
    int indexes_num = linear_idx_norms_mu_U_before_loop_find(I,I2,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,lidx,threshold*residual_pat_norms[0]-residual_pat_norms[1]);


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
 if(cur_frac<0.5){
            indexes_num = lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1]);          
        }
        else{

        indexes_num = lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1]);
        }
        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        t++;

    }
    
    *iters_num_ptr=t-1;
    /*
    printf("threshold=%f frac_for_direct=%f t=%d cur_frac=%f indexes_num=%d\n",threshold,frac_for_direct,t,cur_frac,indexes_num);
    */
    lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);
    return indexes_num;
}





int cauchy_swarz_ncc_linear_idx_clean3(const REAL_TYPE* I,const REAL_TYPE* I2, const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p,REAL_TYPE* norms_of_wins_minus_mu,REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* lidx, REAL_TYPE* vals, REAL_TYPE* U,int best){

    /*struct timeval  tv1, tv2;*/
      
      
printf("XXX cauchy_swarz_ncc_linear_idx_clean3\n");
    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    /*gettimeofday(&tv1,NULL);*/
     
    int indexes_num = linear_idx_norms_mu_U_before_loop_find(I,I2,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,lidx,threshold*residual_pat_norms[0]-residual_pat_norms[1]);

    /*     gettimeofday(&tv2,NULL);*/
    /*double iter_run_time= (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);*/
    /*    printf("%f\n",iter_run_time);*/
    
    box_arr=box_arr+4;
    int total_indexes_num=c_r*c_c;
    REAL_TYPE cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
    int t=1;
    while(t<box_num && cur_frac>frac_for_direct){
        /*  gettimeofday(&tv1,NULL);*/
        int f_r=(*box_arr)-1;
        int t_r=*(box_arr+1)-1;
        int f_c=*(box_arr+2)-1;
        int t_c=*(box_arr+3)-1;
        box_arr=box_arr+4;

        if(cur_frac<0.5){
            indexes_num = lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1]);          
        }
        else{

            indexes_num =  lidx_ncc_inside_loop_scan_all_img(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1],n_p);          
        }
        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        t++;

        /*gettimeofday(&tv2,NULL);*/
        /*        iter_run_time= (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);*/
                  /*printf("%f\n",iter_run_time);*/
   
    }

    /*
    printf("threshold=%f frac_for_direct=%f t=%d cur_frac=%f indexes_num=%d\n",threshold,frac_for_direct,t,cur_frac,indexes_num);
    */
    if(best<0 || indexes_num<best){
        /*
        gettimeofday(&tv1,NULL);
        */
        lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);
        /*
        gettimeofday(&tv2,NULL);
        iter_run_time= (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
        printf("%f\n",iter_run_time);
        */
        return indexes_num;
        
    }
    else{
        /*gettimeofday(&tv1,NULL);*/
        /*printf("best=%d\n",best);*/
        REAL_TYPE mx = U[*lidx];
        int* mx_idx = lidx;
        int* lidx_end=lidx+indexes_num;
        int* lidx_ptr;
        for(lidx_ptr=lidx;lidx_ptr<lidx_end;lidx_ptr++){
            if(U[*lidx_ptr]>mx){
                mx = U[*lidx_ptr];
                mx_idx=lidx_ptr;
            }
        }
        *lidx=*mx_idx;
        lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,1,residual_pat_norms[0],norms_of_wins_minus_mu);

        /*
        gettimeofday(&tv2,NULL);
        iter_run_time= (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
        printf("%f\n",iter_run_time);
        */
        return 1;
        
        /*printf("best=%d\n",best);
        REAL_TYPE mn = U[*lidx];
        int mn_idx = 0;
        int i,j;
        for(i=0;i<best;i++){
            if(U[lidx[i]]<mn){
                mn = U[lidx[i]];
                mn_idx=i;
            }
        }
        for(j=best;j<indexes_num;j++){
            if(U[lidx[j]]>mn){
                lidx[mn_idx]=lidx[j];
                mn = U[*lidx];
                mn_idx=0;
                for(i=0;i<best;i++){
                    if(U[*lidx]<mn){
                        mn = U[*lidx];
                        mn_idx=i;
                    }
                }
                
            }
            }*/
        
        /*
          lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p,
          lidx,vals,best,residual_pat_norms[0],norms_of_wins_minus_mu);
        */
        /*
          return best;
        */
        
    }
    /*    return indexes_num; */
}


int cauchy_swarz_ncc_linear_idx_clean2t(const REAL_TYPE* I,const REAL_TYPE* I2, const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p,REAL_TYPE* norms_of_wins_minus_mu,REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,const REAL_TYPE* thresholds,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* lidx, REAL_TYPE* vals, REAL_TYPE* U,int* iters_num_ptr){

    
    /*    
    int i=0;
    for(i=0;i<20;i++){
        printf("%f ",thresholds[i]);

    }
    */
    /*printf("XXX cauchy_swarz_ncc_linear_idx_clean2t\n");*/
    /*
    for(i=0;i<20;i++){
        printf("%f ",I[i]);

    }
    printf("XXX\n");
    for(i=0;i<20;i++){
        printf("%f ",img[i]);

    }
    
        printf("bu\n");
    */
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;   
    
    int indexes_num =
    linear_idx_norms_mu_U_before_loop_find(I,I2,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p,
    wins_mu,lidx,thresholds[0]*residual_pat_norms[0]-residual_pat_norms[1]);
    /*
    int indexes_num = linear_idx_norms_mu_U_before_loop_find(I,I2,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,lidx,thresholds[0]);
    */
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
        /*        printf("%f\n",thresholds[t]);*/
        if(cur_frac<0.5){
            
            indexes_num =
            lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,thresholds[t]*residual_pat_norms[0]-residual_pat_norms[t+1]);
            /*
            indexes_num =
            lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,thresholds[t]);
            */
        }
        else{
            
            indexes_num =
            lidx_ncc_inside_loop_scan_all_img(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,thresholds[t]*residual_pat_norms[0]-residual_pat_norms[t+1],n_p);
            /*
            indexes_num =  lidx_ncc_inside_loop_scan_all_img(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,thresholds[t],n_p);
            */
            indexes_num =
            lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,thresholds[t]*residual_pat_norms[0]-residual_pat_norms[t+1]);
            
        }
       
        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        t++;
    }
       
    *iters_num_ptr=t-1;
    /*
    printf("threshold=%f frac_for_direct=%f t=%d cur_frac=%f indexes_num=%d\n",threshold,frac_for_direct,t,cur_frac,indexes_num);
    */
    lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);
    return indexes_num;
}


int cauchy_swarz_ncc_linear_idx_clean2t_tl(const REAL_TYPE* I,const REAL_TYPE* I2, const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p,REAL_TYPE* norms_of_wins_minus_mu,REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,const REAL_TYPE* thresholds,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* lidx, REAL_TYPE* vals, REAL_TYPE* U,int* iters_num_ptr,int tl){

    if(tl>0){
        printf("(%f) ",U[tl]);
    }
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;   
    int indexes_num = linear_idx_norms_mu_U_before_loop_find(I,I2,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,lidx,thresholds[0]*residual_pat_norms[0]-residual_pat_norms[1]);

    /*    printf("tl=%d ",tl);*/
    if(tl>0){
        REAL_TYPE u =(U[tl]+residual_pat_norms[1])/residual_pat_norms[0];
        printf("(%f %f) ",u,thresholds[0]);
        printf("(w %f  n %f) ",wins_mu[tl],norms_of_wins_minus_mu[tl]);
        /*        tmp = w*(tmp-support_sz*(*wins_mu))/(*norms_of_wins_minus_mu);*/
        int from_r=box_arr[0]-1;
        int to_r=box_arr[1]-1;
        int from_c=box_arr[2]-1;
        int to_c=box_arr[3]-1;
        int ssz=    (to_r-from_r+1)*(to_c-from_c+1);
        printf("( tmp*w=%f) ",U[tl]*(norms_of_wins_minus_mu[tl])+wins_mu[tl]*ssz*w_arr[0]);
        printf("( w=%f) ",w_arr[0]);
        
        /*        printf("(%f %f) ",U[tl],thresholds[0]*residual_pat_norms[0]-residual_pat_norms[1]);*/
    }
    
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
        if(cur_frac<0.5){
            indexes_num = lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,thresholds[t]*residual_pat_norms[0]-residual_pat_norms[t+1]);          
        }
        else{

            indexes_num = lidx_ncc_inside_loop(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,thresholds[t]*residual_pat_norms[0]-residual_pat_norms[t+1]);
        }
        if(tl>0){
            REAL_TYPE u =(U[tl]+residual_pat_norms[t+1])/residual_pat_norms[0];
            printf("(%f %f) ",u,thresholds[t]);
        }

        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        t++;
    }
    
    if(tl>0){
        printf("\n");
    }
    
    *iters_num_ptr=t-1;
    /*
    printf("threshold=%f frac_for_direct=%f t=%d cur_frac=%f indexes_num=%d\n",threshold,frac_for_direct,t,cur_frac,indexes_num);
    */
    lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);
    return indexes_num;
}



