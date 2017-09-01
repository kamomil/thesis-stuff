#include "cauchy_swarz_ncc_match.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define SEC_RES 1000

int cauchy_swarz_ncc_linear_idx_clean3(const REAL_TYPE* I,const REAL_TYPE* I2, const REAL_TYPE* img,int img_r,int img_c,const REAL_TYPE* pat,int n_p,int m_p,REAL_TYPE* norms_of_wins_minus_mu,REAL_TYPE* wins_mu,const int* box_arr,const REAL_TYPE* w_arr,int box_num,REAL_TYPE threshold,const REAL_TYPE* residual_pat_norms, REAL_TYPE frac_for_direct,int* lidx, REAL_TYPE* vals, REAL_TYPE* U,int best){   
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
            indexes_num =  lidx_ncc_inside_loop_scan_all_img(I,img_r+1,U,norms_of_wins_minus_mu,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],wins_mu,lidx,indexes_num,threshold*residual_pat_norms[0]-residual_pat_norms[t+1],n_p);          
        }
        cur_frac = ((REAL_TYPE) indexes_num)/((REAL_TYPE) total_indexes_num);
        t++;   
    }
    if(best<0 || indexes_num<best){
        lidx_direct_ncc(img,img_r,img_c,pat,n_p,m_p, lidx,vals,indexes_num,residual_pat_norms[0],norms_of_wins_minus_mu);
        return indexes_num;
    }
    else{
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

