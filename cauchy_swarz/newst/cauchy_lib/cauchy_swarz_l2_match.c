#include "cauchy_swarz_l2_match.h"
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>



/*
int cauchy_swarz_inside_loop_opt( img_t_data ,img.rows, img.cols,box_arr_data, box_arr.rows, w_arr_data,threshold,pat_t_data,sz,sz, percent_for_direct, i_loc,j_loc,REAL_TYPE* vals, L);
*/

/*
void cauchy_swarz_inside_loop_opt(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE* L){
*/
    
int cauchy_swarz_inside_loop_opt(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num,const REAL_TYPE* w_arr,REAL_TYPE threshold,const REAL_TYPE* reconst_pat,int n_p,int m_p, REAL_TYPE frac_for_direct,int* i_loc,int* j_loc, REAL_TYPE* vals, REAL_TYPE* L){

    struct timeval  before_loop_tv1, before_loop_tv2;
    struct timeval  loop_tv1, loop_tv2;
    struct timeval  dconv_tv1, dconv_tv2;

    /*printf("CAUCHY_SWARZ_INSIDE_LOOP_OPT\n");
    printf("img: (first 10 ) %f %f %f  %f %f %f  %f %f %f %f\n",img[0],img[1],img[2],img[3],img[4],img[5],img[6],img[7],img[8],img[9]);
    printf("img_r=%d, img_c=%d\n",img_r,img_c);
    printf("box_arr first 5 %d %d %d %d %d\n",box_arr[0],box_arr[1],box_arr[2],box_arr[3],box_arr[4]);
    printf("box_num=%d, threshold=%f, n_p=%d, m_p=%d\n",box_num,threshold,n_p,m_p);
    printf("w_arr (first 5) %f %f %f %f %f\n",w_arr[0],w_arr[1],w_arr[2],w_arr[3],w_arr[4]);
    printf("pat: (first 3 ) %f %f %f \n",reconst_pat[0],reconst_pat[1],reconst_pat[2]);*/
                    
    gettimeofday(&before_loop_tv1,NULL);
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

#ifdef MEX
    REAL_TYPE* I =  mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
    REAL_TYPE* I_sqrs =  mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
    REAL_TYPE* wins_ss_plus_pat_ss= mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* wins_norms= mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* C= mxMalloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* pat_cpy= mxMalloc(sizeof(REAL_TYPE)*n_p*m_p);
    
#else
    REAL_TYPE* I = /*(REAL_TYPE*)*/ malloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
    REAL_TYPE* I_sqrs = /*(REAL_TYPE*)*/ malloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
    REAL_TYPE* wins_ss_plus_pat_ss = /*(REAL_TYPE*)*/ malloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* wins_norms = /*(REAL_TYPE*)*/ malloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* C = /*(REAL_TYPE*)*/ malloc(sizeof(REAL_TYPE)*c_r*c_c);
    REAL_TYPE* pat_cpy = malloc(sizeof(REAL_TYPE)*n_p*m_p);
    
    /*    REAL_TYPE* I = new REAL_TYPE[sizeof(REAL_TYPE)*(img_r+1)*(img_c+1)];
    REAL_TYPE* I_sqrs = new REAL_TYPE[sizeof(REAL_TYPE)*(img_r+1)*(img_c+1)];
    REAL_TYPE* wins_sum_sqrs = new REAL_TYPE[sizeof(REAL_TYPE)*c_r*c_c];
    REAL_TYPE* wins_norms = new REAL_TYPE[sizeof(REAL_TYPE)*c_r*c_c];
    REAL_TYPE* C = new REAL_TYPE[sizeof(REAL_TYPE)*c_r*c_c];*/
    
#endif

    memcpy(pat_cpy,reconst_pat,sizeof(REAL_TYPE)*n_p*m_p);
    int_img_and_int_img_sqr(img, img_r, img_c ,I,  I_sqrs);
        

    REAL_TYPE pat_sum_sqrs =sum_sqrs(reconst_pat,n_p,m_p);
    

    /*    L_before_loop(I_sqrs,img_r,img_c,n_p,m_p, mul2,pat_sum_sqrs, wins_sum_sqrs, wins_norms, L);*/
    /*    free(I_sqrs);*/
    
    threshold = threshold*n_p*m_p;
        
    int t=0;

    int* i[2];
    int* j[2];

#ifdef MEX

    i[0] = mxMalloc(c_r*c_c*sizeof(int));
    j[0] = mxMalloc(c_r*c_c*sizeof(int));
    i[1] = mxMalloc(c_r*c_c*sizeof(int));
    j[1] = mxMalloc(c_r*c_c*sizeof(int));
    
#else
    i[0]=/*(int*)*/malloc(c_r*c_c*sizeof(int));
    j[0]=/*(int*)*/malloc(c_r*c_c*sizeof(int));
    i[1]=/*(int*)*/malloc(c_r*c_c*sizeof(int));
    j[1]=/*(int*)*/malloc(c_r*c_c*sizeof(int));

    /*i[0]=new int[c_r*c_c*sizeof(int)]
    j[0]=new int[c_r*c_c*sizeof(int)];
    i[1]=new int[c_r*c_c*sizeof(int)];
    j[1]=new int[c_r*c_c*sizeof(int)];*/
    
#endif
    int indicator=0;
    /*    int indices_num = find2(L, c_r,c_c, threshold,i[indicator],j[indicator]);*/
    int indices_num = L_before_loop_and_find(I_sqrs,img_r,img_c,n_p,m_p,pat_sum_sqrs, wins_ss_plus_pat_ss, wins_norms, L,i[indicator],j[indicator],threshold);    
    int total_indices=c_r*c_c;
    REAL_TYPE cur_frac = ((REAL_TYPE) indices_num)/((REAL_TYPE) total_indices);
        
#ifdef MEX
    mxFree(I_sqrs);
#else
    free(I_sqrs);
        /*delete[] I_sqrs;    */
#endif
    
    zero2(C,c_c*c_r);
    REAL_TYPE residaul_pat_ss =pat_sum_sqrs;


    
    gettimeofday(&before_loop_tv2,NULL);
    gettimeofday(&loop_tv1,NULL);
    while(t<box_num && cur_frac>frac_for_direct*10){//TODO - remove *10
        /*printf("frc is %f\n",cur_frac);*/
        int f_r=(*box_arr)-1;
        int t_r=*(box_arr+1)-1;
        int f_c=*(box_arr+2)-1;
        int t_c=*(box_arr+3)-1;
        box_arr=box_arr+4;
        int support_sz=(t_r-f_r+1)*(t_c-f_c+1);

        //residual_p_ss at the t'th iteration is the ss of
        //(P-sum(box(1)..box(t))), (TODO - check that it works)
        residaul_pat_ss=residaul_pat_ss + (w_arr[t]*w_arr[t]*support_sz) - (2*w_arr[t]*remove_box2(pat_cpy,n_p,m_p,f_r,t_r,f_c,t_c,w_arr[t])) ;
        
        /*printf("residual pat sumsqrs       =%f\n", residaul_pat_ss);
          printf("reconst pat sumsqrs after  =%f\n\n", sum_sqrs(reconst_pat,n_p,m_p));*/
        indices_num = inside_loop3_tmp(I,img_r+1,C,L,wins_ss_plus_pat_ss,wins_norms,c_r,c_c,f_r,t_r,f_c,t_c,w_arr[t],pat_sum_sqrs,i[indicator], j[indicator],i[1-indicator],j[1-indicator], indices_num,threshold, residaul_pat_ss);
        t++;
        indicator=1-indicator;
        cur_frac = ((REAL_TYPE) indices_num)/((REAL_TYPE) total_indices);
    }


    
    gettimeofday(&loop_tv2,NULL);

    int* end=i[1-indicator]+indices_num;
    int* i_ptr=i[1-indicator];
    int* j_ptr=j[1-indicator];
    
    int* i_loc_ptr=i_loc;
    int* j_loc_ptr=j_loc;
    
    for (;i_ptr<end;i_ptr++,j_ptr++,i_loc_ptr++,j_loc_ptr++){
        *i_loc_ptr = *i_ptr;
        *j_loc_ptr = *j_ptr;
    }

    /*
   
    printf("indices_num=%d\n",indices_num);
    printf("i[1-indicator] %d %d %d\n",i[1-indicator][0], i[1-indicator][indices_num-1],i[1-indicator][300]);
    printf("j[1-indicator] %d %d %d\n",j[1-indicator][0], j[1-indicator][indices_num-1],j[1-indicator][300]);
    printf("i_loc %d %d %d\n",i_loc[0], i_loc[indices_num-1],i_loc[300]);
    printf("j_loc %d %d %d\n",j_loc[0], j_loc[indices_num-1],j_loc[300]);
    */

    //TODO - this is just for experiments for measuring time for different
    //indexes percent

    int num_to_remove= (int) indices_num-frac_for_direct*total_indices+2;
    indices_num=indices_num-num_to_remove;
    cur_frac = ((REAL_TYPE) indices_num)/((REAL_TYPE) total_indices);
    
    //////////////////////////////////////
    //////////////////////////////////////
    
    gettimeofday(&dconv_tv1,NULL);
    direct_conv(img,img_r,img_c,reconst_pat,n_p,m_p, i_loc, j_loc,vals,indices_num);
    gettimeofday(&dconv_tv2,NULL);
    //printf("vals %f %f %f\n",vals[0], vals[indices_num-1],vals[300]);
    
#ifdef MEX
    mxFree(I);
    mxFree(wins_sum_sqrs);
    mxFree(wins_norms);
    mxFree(C);
    mxFree(i[0]);
    mxFree(j[0]);
    mxFree(i[1]);
    mxFree(j[1]);
    mxFree(pat_cpy);

#else
    free(I);
    free(wins_ss_plus_pat_ss);
    free(wins_norms);
    free(C);
    free(i[0]);
    free(j[0]);
    free(i[1]);
    free(j[1]);
    free(pat_cpy);
    
    /*delete[] I;
    delete[] wins_sum_sqrs;
    delete[] wins_norms;
    delete[] C;
    delete[] i[0];
    delete[] j[0];
    delete[] i[1];
    delete[] j[1];*/
#endif

    //printf ("Total Box time before loop = %f seconds, bunm=%d\n",(double) (before_loop_tv2.tv_usec - before_loop_tv1.tv_usec)/1000000 + (double) (before_loop_tv2.tv_sec - before_loop_tv1.tv_sec),box_num);
    //printf ("Total Box time inside loop = %f seconds, bunm=%d\n",(double) (loop_tv2.tv_usec - loop_tv1.tv_usec) / 1000000 +(double) (loop_tv2.tv_sec - loop_tv1.tv_sec),box_num);
    printf ("Total Box time direct conv = %f seconds, frac=%f\n",(double) (dconv_tv2.tv_usec - dconv_tv1.tv_usec) / 1000000 +(double) (dconv_tv2.tv_sec - dconv_tv1.tv_sec),cur_frac);

    return indices_num;
}
