#include "box_corr_lib.h"


#ifdef MEX 
#include "matrix.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void indices_mat2c(int* i,int size){

    int k;
    for(k=0;k<size;k++){
        i[k]--;
    }
    
}
void integral_image_opt(const REAL_TYPE* img,int img_r, int img_c ,REAL_TYPE* I){

    
    int i=0;
    REAL_TYPE* I_end_ptr = I+(img_r+1)*(img_c+1);
    REAL_TYPE* I_end_col_ptr=I+img_r+1;
    
    REAL_TYPE col_sum=0;
    REAL_TYPE* I_ptr=I;
    REAL_TYPE* I_ptr2=I;
    const REAL_TYPE* img_ptr=img;
    
    for(I_ptr=I;I_ptr<I_end_col_ptr;I_ptr++){/*zero the first column*/
        *I_ptr=0;    
    }

    while(I_ptr<I_end_ptr){
        *I_ptr=0;
        col_sum=0;
        for(i=0;i<img_r;i++){
            I_ptr++;
            I_ptr2++;
            col_sum=col_sum+(*img_ptr);
            *I_ptr=(*I_ptr2)+col_sum;
            img_ptr++;
        }
        I_ptr++;
        I_ptr2++;
    }
}

void int_img_and_int_img_sqr(const REAL_TYPE* img, int img_r, int img_c ,REAL_TYPE* I, REAL_TYPE* I2){
    

    /*printf("int_img_and_int_img_sqr\n");*/
    int i=0;
    int I_size = (img_r+1)*(img_c+1);
    int I_r=img_r+1;
    
    REAL_TYPE col_sum_I=0;
    REAL_TYPE col_sum_I2=0;
    REAL_TYPE* I_ptr=I;
    REAL_TYPE* I_ptr2=I;
    REAL_TYPE* I2_ptr=I2;
    REAL_TYPE* I2_ptr2=I2;
    
    const REAL_TYPE* img_ptr=img;
    
    for(;I_ptr<I+I_r;I_ptr++,I2_ptr++){/*zero the first column*/
        *I_ptr=0;
        *I2_ptr=0;    
    }

    while(I_ptr<I+I_size){
        *I_ptr=0;
        *I2_ptr=0;
        col_sum_I=0;
        col_sum_I2=0;
        for(i=0;i<img_r;i++){
            I_ptr++;
            I_ptr2++;
            col_sum_I=col_sum_I+(*img_ptr);
            *I_ptr=(*I_ptr2)+col_sum_I;

            I2_ptr++;
            I2_ptr2++;
            col_sum_I2=col_sum_I2+((*img_ptr)*(*img_ptr));
            *I2_ptr=(*I2_ptr2)+col_sum_I2;
            
            img_ptr++;
        }
        I_ptr++;
        I_ptr2++;

        I2_ptr++;
        I2_ptr2++;
    }
}


/*
##########################
#####
#########################
 */


int inside_loop(const REAL_TYPE* I,int I_r,REAL_TYPE* C, REAL_TYPE* L,const REAL_TYPE* wins_sum_sqrs,const REAL_TYPE* wins_norms,int c_r,int c_c,int from_r ,int to_r ,int from_c ,int to_c, REAL_TYPE w, REAL_TYPE p_ss, REAL_TYPE* reconst_pat,int n_p,int m_p, const int* i_cur, const int* j_cur, int* i_next, int* j_next,int ind_num,REAL_TYPE threshold){

    int num=0;
    int idxidx=0;
    int ci,Ij,Ii,cj=0;

    
    remove_box(reconst_pat,n_p,m_p,from_r,to_r,from_c,to_c,w);
    REAL_TYPE mul=-2*sqrt(sum_sqrs(reconst_pat,n_p,m_p));

    while(cj<c_c && idxidx<ind_num){
        cj=j_cur[idxidx];
        Ij=cj+1;
        
        while(j_cur[idxidx]==cj && idxidx<ind_num){
            ci=i_cur[idxidx];
            Ii=ci+1;
            
            int idx = ci+c_r*cj;
            
            REAL_TYPE tmp = I[(Ii+to_r)+I_r*(Ij+to_c)];
            tmp=tmp-I[(Ii+to_r)+I_r*(Ij+from_c-1)];
            tmp=tmp-I[(Ii+from_r-1)+I_r*(Ij+to_c)];
            tmp=tmp+I[(Ii+from_r-1)+I_r*(Ij+from_c-1)];

            C[idx]=C[idx]+w*tmp;
            L[idx]=1.0*wins_sum_sqrs[idx]+(mul)*wins_norms[idx]+(-2)*C[idx]+p_ss;
            if(L[idx]<threshold){
                i_next[num]=i_cur[idxidx];
                j_next[num]=j_cur[idxidx];
                num++;
            }
            idxidx++;
        }       
    }
    return num;
}  


int inside_loop3(const REAL_TYPE* I,int I_r,REAL_TYPE* C, REAL_TYPE* L,const REAL_TYPE* wins_ss_plus_pat_ss,const REAL_TYPE* wins_norms,int c_r,int c_c,int from_r ,int to_r ,int from_c ,int to_c, REAL_TYPE w, REAL_TYPE p_ss, REAL_TYPE* reconst_pat,int n_p,int m_p, const int* i_cur_ptr, const int* j_cur_ptr, int* i_next_ptr, int* j_next_ptr,int ind_num,REAL_TYPE threshold){

    int num=0;
    int ci,cj;

    
    /*TODO, this can be more efficient, by calculating the new pat sum squers
      from the previous one */
    remove_box2(reconst_pat,n_p,m_p,from_r,to_r,from_c,to_c,w);

    REAL_TYPE mul=-2*sqrt(sum_sqrs(reconst_pat,n_p,m_p));
    const int* end_i_cur_ptr=i_cur_ptr+ind_num;
    
    const REAL_TYPE* I_br_ptr;
    const REAL_TYPE* I_bl_ptr;
    const REAL_TYPE* I_ur_ptr;
    const REAL_TYPE* I_ul_ptr;
 
    const REAL_TYPE* ul_const_term = I+from_r + I_r*from_c;
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);

    REAL_TYPE* C_ptr;
    
    while(i_cur_ptr<end_i_cur_ptr){
        cj=(*j_cur_ptr);

       
        int I_col_term = I_r*cj;
        int C_col_term = c_r*cj;
        
        
        while( (*j_cur_ptr)==cj && i_cur_ptr<end_i_cur_ptr){
            ci=(*i_cur_ptr);
              
            I_ul_ptr = ul_const_term+ci+I_col_term;
            I_ur_ptr = ur_const_term+ci+I_col_term;
            I_bl_ptr = bl_const_term+ci+I_col_term;
            I_br_ptr = br_const_term+ci+I_col_term;
            
            REAL_TYPE tmp = (*I_br_ptr);
            tmp -=(*I_ur_ptr);
            tmp -=(*I_bl_ptr);
            tmp +=(*I_ul_ptr);

            C_ptr=C+ci+C_col_term;
            (*C_ptr) += (w*tmp);
            (*(L+ci+C_col_term))=(*(wins_ss_plus_pat_ss+ci+C_col_term))+(mul)*(*(wins_norms+ci+C_col_term))-2*(*C_ptr);
            if(*(L+ci+C_col_term)<threshold){
                (*i_next_ptr)=(*i_cur_ptr);
                (*j_next_ptr)=(*j_cur_ptr);
                i_next_ptr++;
                j_next_ptr++;
                num++;
            }
            i_cur_ptr++;
            j_cur_ptr++;
        }
        
    }
    return num;
}  



int inside_loop3_tmp(const REAL_TYPE* I,int I_r,REAL_TYPE* C, REAL_TYPE* L,const REAL_TYPE* wins_ss_plus_pat_ss,const REAL_TYPE* wins_norms,int c_r,int c_c,int from_r ,int to_r ,int from_c ,int to_c, REAL_TYPE w, REAL_TYPE p_ss,  const int* i_cur_ptr, const int* j_cur_ptr, int* i_next_ptr, int* j_next_ptr,int ind_num,REAL_TYPE threshold,REAL_TYPE residual_p_ss){

    int num=0;
    int ci,cj;

    REAL_TYPE mul=-2*sqrt(residual_p_ss);
    const int* end_i_cur_ptr=i_cur_ptr+ind_num;    
    const REAL_TYPE* I_br_ptr;
    const REAL_TYPE* I_bl_ptr;
    const REAL_TYPE* I_ur_ptr;
    const REAL_TYPE* I_ul_ptr;
    const REAL_TYPE* ul_const_term = I+from_r + I_r*from_c;
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);

    REAL_TYPE* C_ptr;
    REAL_TYPE tmp;

    int I_col_term;
    int C_col_term;
        
    while(i_cur_ptr<end_i_cur_ptr){
        cj=(*j_cur_ptr);

        I_col_term = I_r*cj;
        C_col_term = c_r*cj;
        
        while( (*j_cur_ptr)==cj && i_cur_ptr<end_i_cur_ptr){
            ci=(*i_cur_ptr);
              
            I_ul_ptr = ul_const_term+ci+I_col_term;
            I_ur_ptr = ur_const_term+ci+I_col_term;
            I_bl_ptr = bl_const_term+ci+I_col_term;
            I_br_ptr = br_const_term+ci+I_col_term;
            
            tmp = (*I_br_ptr);
            tmp -=(*I_ur_ptr);
            tmp -=(*I_bl_ptr);
            tmp +=(*I_ul_ptr);

            C_ptr=C+ci+C_col_term;
            (*C_ptr) += (w*tmp);
            REAL_TYPE tmp2 = (*(wins_ss_plus_pat_ss+ci+C_col_term))+(mul)*(*(wins_norms+ci+C_col_term))-2*(*C_ptr);
            /*(*(L+ci+C_col_term))=(*(wins_ss_plus_pat_ss+ci+C_col_term))+(mul)*(*(wins_norms+ci+C_col_term))-2*(*C_ptr);*/
            if(tmp2<threshold){
                (*(L+ci+C_col_term))=tmp2;
                /*if(*(L+ci+C_col_term)<threshold){*/
                (*i_next_ptr)=(*i_cur_ptr);
                (*j_next_ptr)=(*j_cur_ptr);
                i_next_ptr++;
                j_next_ptr++;
                num++;
            }
            i_cur_ptr++;
            j_cur_ptr++;
        }
    }
    return num;
}  

/*
##########################
#####
#########################
*/

void box_corr_matlab_indices(const REAL_TYPE* img,int img_r,int img_c, const int* box_arr,int box_num, const REAL_TYPE* w_arr,int n_p,int m_p, REAL_TYPE* C){    

#ifdef MEX
    REAL_TYPE* I = mxMalloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
#else
    REAL_TYPE* I = /*(REAL_TYPE*)*/ malloc(sizeof(REAL_TYPE)*(img_r+1)*(img_c+1));
#endif
    
    integral_image_opt(img,img_r,img_c,I);

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    int bidx=0;

    int ci,cj,Ii,Ij;

    zero(C,c_r*c_c);
 

    for(bidx=0;bidx<box_num;bidx++){
 
        int from_r = box_arr[bidx]-1;
        int to_r = box_arr[bidx+box_num]-1;
        int from_c = box_arr[bidx+2*box_num]-1;
        int to_c = box_arr[bidx+3*box_num]-1;
   
        Ii=1;
        for(ci=0;ci<c_c;ci++){
            Ij=1;
            for(cj=0;cj<c_r;cj++){

                /*                printf("br=%d\n",(Ij+to_r)+(img_r+1)*(Ii+to_c));
                printf("ur=%d\n",(Ij+from_r-1)+(img_r+1)*(Ii+to_c));
                printf("bl=%d\n",(Ij+to_r)+(img_r+1)*(Ii+from_c-1));
                printf("ul=%d\n",(Ij+from_r-1)+(img_r+1)*(Ii+from_c-1));*/

                
                REAL_TYPE tmp = I[(Ij+to_r)+(img_r+1)*(Ii+to_c)];
                tmp=tmp-I[(Ij+to_r)+(img_r+1)*(Ii+from_c-1)];/*bl*/
                tmp=tmp-I[(Ij+from_r-1)+(img_r+1)*(Ii+to_c)];
                tmp=tmp+I[(Ij+from_r-1)+(img_r+1)*(Ii+from_c-1)];
                C[c_r*ci+cj]=C[c_r*ci+cj]+w_arr[bidx]*tmp;
                Ij++;
            }
            Ii++;
        }     
    }
#ifdef MEX
    mxFree(I);
#else
    free(I);
#endif
}
   



void wins_sums(const REAL_TYPE* I,int img_r,int img_c, int n_p,int m_p, REAL_TYPE* C){    

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int I_r=img_r+1;
    
    REAL_TYPE* C_ptr=C;
    REAL_TYPE* C_end_ptr=C+c_r*c_c;
    
    const REAL_TYPE* I_br_ptr=I+I_r*m_p+n_p;
    const REAL_TYPE* I_bl_ptr=I+n_p;
    const REAL_TYPE* I_ur_ptr=I+I_r*m_p;
    const REAL_TYPE* I_ul_ptr=I;

    while(C_ptr<C_end_ptr){
        REAL_TYPE* next_col_ptr = C_ptr+c_r;
        while(C_ptr<next_col_ptr){

            REAL_TYPE tmp = (*I_br_ptr);
            tmp -=(*I_ur_ptr);
            tmp -=(*I_bl_ptr);
            tmp +=(*I_ul_ptr);
            (*C_ptr)=tmp;

            I_br_ptr++;
            I_ur_ptr++;
            I_bl_ptr++;
            I_ul_ptr++;
            C_ptr++;            
        }

        I_br_ptr=I_br_ptr+n_p;
        I_ur_ptr=I_ur_ptr+n_p;
        I_bl_ptr=I_bl_ptr+n_p;
        I_ul_ptr=I_ul_ptr+n_p;
    }
}

void L_before_loop(const REAL_TYPE* I,int img_r,int img_c, int n_p,int m_p, REAL_TYPE mul2,REAL_TYPE p_ss, REAL_TYPE* ss_ptr, REAL_TYPE* norms_ptr, REAL_TYPE* L_ptr){    

    /*

    wins_sums(I_sqrs,img_r,img_c, n_p, m_p,wins_sum_sqrs);
    mxFree(I_sqrs);

    matrix_sqrt2(wins_sum_sqrs,c_r*c_c, wins_norms);

    void matrix_add(const REAL_TYPE* m1, const REAL_TYPE* m2,REAL_TYPE addScalar, REAL_TYPE mulScalar1, REAL_TYPE mulScalar2, int size, REAL_TYPE* out){
    matrix_add2(wins_sum_sqrs, wins_norms, pat_sum_sqrs , 1.0,mul2, c_r*c_c, L);

     */
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int I_r=img_r+1;
    
    REAL_TYPE* ss_end_ptr=ss_ptr+c_r*c_c;
    
    const REAL_TYPE* I_br_ptr=I+I_r*m_p+n_p;
    const REAL_TYPE* I_bl_ptr=I+n_p;
    const REAL_TYPE* I_ur_ptr=I+I_r*m_p;
    const REAL_TYPE* I_ul_ptr=I;

    while(ss_ptr<ss_end_ptr){
        REAL_TYPE* next_col_ptr = ss_ptr+c_r;
        while(ss_ptr<next_col_ptr){

            REAL_TYPE tmp = (*I_br_ptr);
            tmp -=(*I_ur_ptr);
            tmp -=(*I_bl_ptr);
            tmp +=(*I_ul_ptr);

            REAL_TYPE sq_tmp=sqrt(tmp);
            
            (*ss_ptr)=tmp;
            (*norms_ptr)=sq_tmp;
            (*L_ptr)=tmp+mul2*sq_tmp+p_ss;
                
            I_br_ptr++;
            I_ur_ptr++;
            I_bl_ptr++;
            I_ul_ptr++;
            ss_ptr++;
            norms_ptr++;
            L_ptr++;
        }

        I_br_ptr=I_br_ptr+n_p;
        I_ur_ptr=I_ur_ptr+n_p;
        I_bl_ptr=I_bl_ptr+n_p;
        I_ul_ptr=I_ul_ptr+n_p;
    }
}

/*    int indices_num = L_before_loop_and_find(I_sqrs,img_r,img_c,n_p,m_p, mul2,pat_sum_sqrs, wins_sum_sqrs, wins_norms, L,i[indicator],j[indicator],threshold);    */
int L_before_loop_and_find(const REAL_TYPE* I_sqrs,int img_r,int img_c, int n_p,int m_p, REAL_TYPE p_ss, REAL_TYPE* ss_ptr, REAL_TYPE* norms_ptr, REAL_TYPE* L_ptr,int* i_ptr,int* j_ptr,REAL_TYPE threshold){

    REAL_TYPE mul2=2*sqrt(p_ss);

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int I_r=img_r+1;
    int row=0;
    int col=0;
    int num=0;
    REAL_TYPE* ss_end_ptr=ss_ptr+c_r*c_c;
    
    const REAL_TYPE* I_br_ptr=I_sqrs+I_r*m_p+n_p;
    const REAL_TYPE* I_bl_ptr=I_sqrs+n_p;
    const REAL_TYPE* I_ur_ptr=I_sqrs+I_r*m_p;
    const REAL_TYPE* I_ul_ptr=I_sqrs;

    while(ss_ptr<ss_end_ptr){
        REAL_TYPE* next_col_ptr = ss_ptr+c_r;
        while(ss_ptr<next_col_ptr){

            REAL_TYPE tmp = (*I_br_ptr);
            tmp -=(*I_ur_ptr);
            tmp -=(*I_bl_ptr);
            tmp +=(*I_ul_ptr);

            REAL_TYPE sq_tmp=sqrt(tmp);
            
            (*ss_ptr)=tmp+p_ss;
            (*norms_ptr)=sq_tmp;
            (*L_ptr)=tmp-mul2*sq_tmp+p_ss;
            if(*L_ptr<threshold){
                *i_ptr=row;
                *j_ptr=col;
                i_ptr++;
                j_ptr++;
                num++;
            }           
            I_br_ptr++;
            I_ur_ptr++;
            I_bl_ptr++;
            I_ul_ptr++;
            ss_ptr++;
            norms_ptr++;
            L_ptr++;
            row++;
        }
        col++;
        row=0;
        
        I_br_ptr=I_br_ptr+n_p;
        I_ur_ptr=I_ur_ptr+n_p;
        I_bl_ptr=I_bl_ptr+n_p;
        I_ul_ptr=I_ul_ptr+n_p;
    }
    return num;
}


/*
######################3
####
#######################
 */

REAL_TYPE sum_sqrs(const REAL_TYPE* p,int n_p,int m_p){

    REAL_TYPE s=0;
    int i;
    for(i=0;i<n_p*m_p;i++){
        s=s+(p[i]*p[i]);
    }
    /*      printf("IN SUM_SQRS %f ",s);*/
    return s;
}

void matrix_affine_trans(const REAL_TYPE* m, REAL_TYPE addScalar, REAL_TYPE mulScalar, int size, REAL_TYPE* out){

    
    int i;
    for(i=0;i<size;i++){
        /*        printf("mat aff i=%d s=%d",i,size);*/
        out[i]= mulScalar*m[i]+addScalar;
    }
}

void matrix_add(const REAL_TYPE* m1, const REAL_TYPE* m2,REAL_TYPE addScalar, REAL_TYPE mulScalar1, REAL_TYPE mulScalar2, int size, REAL_TYPE* out){

    int i;
    for(i=0;i<size;i++){        
        out[i]= mulScalar1*m1[i]+mulScalar2*m2[i]+addScalar;
    }
}

void matrix_add2(const REAL_TYPE* m1, const REAL_TYPE* m2,REAL_TYPE addScalar, REAL_TYPE mulScalar1, REAL_TYPE mulScalar2, int size, REAL_TYPE* out){

    REAL_TYPE* to=out+size;
    for(;out<to;out++,m1++,m2++){        
        (*out) = mulScalar1*(*m1) + mulScalar2*(*m2) + addScalar;
    }
}

void matrix_sqrt(const REAL_TYPE* m, int size, REAL_TYPE* out){

    int i;
    for(i=0;i<size;i++){
        out[i]= sqrt(m[i]);
    }
}

void matrix_sqrt2(const REAL_TYPE* m, int size, REAL_TYPE* out){

    REAL_TYPE* to = out+size;
    for(;out<to;out++,m++){
        (*out)= sqrt(*m);
    }
}

int find_from_indices(const REAL_TYPE* m,int dim1,int dim2,REAL_TYPE upperbound,int ind_num, int * i_cur,int* j_cur, int* i_next, int* j_next){

    int k;
    int num=0;
    for(k=0;k<ind_num;k++){
        if(m[i_cur[k]+dim1*j_cur[k]]<upperbound){
            i_next[num]=i_cur[k];
            j_next[num]=j_cur[k];
            num++;
        }
    }
    return num;
}

int find(const REAL_TYPE* m, int dim1, int dim2, REAL_TYPE upperbound, int* i, int* j){

    int idx,jdx,jcount,num=0;
    jcount=0;
    for(jdx=0;jdx<dim1*dim2;jdx+=dim1){
        for (idx=0;idx<dim1;idx++){
            if(m[idx+jdx]<upperbound){
                i[num]=idx;
                j[num]=jcount;
                num++;
            }
        }
        jcount++;
    }
    return num;
}


int find2(const REAL_TYPE* m, int dim1, int dim2, REAL_TYPE upperbound, int* i, int* j){

    const REAL_TYPE* end_ptr=m+dim1*dim2;

    int row=0;
    int col=0;
    int num=0;
    
    while(m<end_ptr){
        const REAL_TYPE* m_next_col=m+dim1;
        for(;m<m_next_col;m++){
            if(*m<upperbound){
                *i=row;
                *j=col;
                i++;
                j++;
                num++;
            }
            row++;
        }
        col++;
        row=0;
    }
    return num;
}


int find_new(REAL_TYPE p_norm,const REAL_TYPE* wins_norms,const REAL_TYPE* thresh_minus_const_term, int c_r, int c_c, int* i,int* j){

    const REAL_TYPE* wins_norms_ptr = wins_norms;
    const REAL_TYPE* thresh_ptr = thresh_minus_const_term;

    int* i_ptr=i;
    int* j_ptr=j;

    int idx=0;
    int jdx=0;
    int num=0;
    for(;thresh_ptr<thresh_minus_const_term+(c_r*c_c);thresh_ptr++,wins_norms_ptr++){
        
        if(-1*(*wins_norms_ptr)*p_norm < *thresh_ptr){
            *i_ptr=idx;
            *j_ptr=jdx;
            num++;
            
        }
    
        idx++;
        if(idx==(jdx+1)*c_r){
            idx=0;
            jdx++;
        }
    }
    
    return num;
}
                                                                                                                               

  /* f=zeros(n_p,m_p);
         * f(box_arr(t,1):box_arr(t,2),box_arr(t,3):box_arr(t,4))=w_arr(t);
         * reconst_filter=reconst_filter-f;
         */
      
void remove_box(REAL_TYPE* pat,int n_p,int m_p, int b1,int b2,int b3,int b4,REAL_TYPE c){

    int idx,jdx;
    for(jdx=n_p*b3;jdx<=n_p*b4;jdx+=n_p){/* run on columns */
        for (idx=b1;idx<=b2;idx++){ /* run on rows */
            pat[idx+jdx]=pat[idx+jdx]-c;
        }
    }
}

/*remove the box from the pattern and return the sum of pattern's elements in the support
  of thebox*/
REAL_TYPE remove_box2(REAL_TYPE* pat,int n_p,int m_p, int f_r,int t_r,int f_c,int t_c,REAL_TYPE c){

    
    REAL_TYPE* pat_ptr=pat+(n_p*f_c+f_r);
    REAL_TYPE* end_ptr=pat+(n_p*t_c+t_r);
    REAL_TYPE sum=0;
    while(pat_ptr<end_ptr){
        REAL_TYPE* end_col_ptr=pat_ptr+t_r-f_r+1;
        REAL_TYPE* start_next_col_ptr=pat_ptr+n_p;
        for(;pat_ptr<end_col_ptr;pat_ptr++){
            sum+=(*pat_ptr);
            (*pat_ptr)=(*pat_ptr)-c;
        }
        pat_ptr=start_next_col_ptr;
    }
    return sum;
}

void remove_box_matlab_indices(REAL_TYPE* pat,int n_p,int m_p, int b1,int b2,int b3,int b4,REAL_TYPE c){

    
    int idx,jdx;
    for(jdx=n_p*(b3-1);jdx<=n_p*(b4-1);jdx+=n_p){/* run on columns */
        for (idx=b1-1;idx<=b2-1;idx++){ /* run on rows */
            pat[idx+jdx]=pat[idx+jdx]-c;
        }
    }
}
   

void zero(REAL_TYPE* m,int size){

    int i;
    for(i=0;i<size;i++){
        m[i]=0;
    }


}

void zero2(REAL_TYPE* m,int size){

    REAL_TYPE* to = m+size;
    
    for(;m<to;m++){
        (*m)=0;
    }


}

void get_vals_in_indices(const REAL_TYPE* m,int dim1,int dim2,const int* i, const int* j,int ind_num,REAL_TYPE* vals){

    int k;
    for(k=0;k<ind_num;k++){
        vals[k]=m[i[k]+dim1*j[k]];

    }
}

void set_vals_in_indices(REAL_TYPE* m,int dim1,const int* i,const int* j,int ind_num,const REAL_TYPE* vals){

    int k;
    for(k=0;k<ind_num;k++){
        m[i[k]+dim1*j[k]]=vals[k];
    }

}


void mat_printf(const REAL_TYPE* m,int dim1,int from_r,int to_r,int from_c,int to_c){

    int idx,jdx;
    for (idx=from_r;idx<=to_r;idx++){ /* run on rows */
        for(jdx=dim1*from_c;jdx<=dim1*to_c;jdx+=dim1){/* run on columns */
            /*printf("%d ",m[idx+jdx]);*/
            printf("%f ",m[idx+jdx]);
        }
        printf("\n");
        
    }
}
/*   add_3_matrices_in_indices(wins_sum_sqrs,1.0,           wins_norms,   mul,C,-2,pat_sum_sqrs,c_r,i,j,indices_num,L);*/
void add_3_matrices_in_indices(const REAL_TYPE* m1,REAL_TYPE mul1,const REAL_TYPE* m2,REAL_TYPE mul2,const REAL_TYPE* m3,REAL_TYPE mul3,REAL_TYPE adder,int dim1,const int* i, const int* j,int ind_num,REAL_TYPE* L){

    int k;
    for(k=0;k<ind_num;k++){
        int idx = i[k]+dim1*j[k];   
        L[idx]=mul1*m1[idx]+mul2*m2[idx]+mul3*m3[idx]+adder;
    }
}

REAL_TYPE mat_dist(const REAL_TYPE* m1, const REAL_TYPE* m2, int s){

    REAL_TYPE sum=0;
    int k;
    for(k=0;k<s;k++){
        sum=sum+(m1[k]-m2[k])*(m1[k]-m2[k]);
    }
    return sum;
}


REAL_TYPE max_diff(const REAL_TYPE* m1, const REAL_TYPE* m2, int s){

    REAL_TYPE max=0;
    int k;
    for(k=0;k<s;k++){
        if(abs(m1[k]-m2[k])>max){
            max=abs(m1[k]-m2[k]);
        }
    }
    return max;
}
    
void mat_copy(const REAL_TYPE* f, REAL_TYPE* t, int s){

    int k;
    for(k=0;k<s;k++){
        t[k]=f[k];
    }
}


void reconstruct_pat(const int* box_arr,int box_num, const REAL_TYPE* w_arr,int n_p,int m_p, REAL_TYPE* pat){

    int i;
    zero2(pat,n_p*m_p);


    /*for(i=0;i<1;i++){*/
    for(i=0;i<box_num;i++){
        int f_r = box_arr[i]-1;
        int t_r = box_arr[i+box_num]-1;
        int f_c = box_arr[i+2*box_num]-1;
        int t_c = box_arr[i+3*box_num]-1;
        
        REAL_TYPE c = w_arr[i];
        /*  printf("%d %d %d %d %f\n",f_r,t_r,f_c,t_c,c);*/
 
        REAL_TYPE* pat_ptr=pat+(n_p*f_c+f_r);
        REAL_TYPE* end_ptr=pat+(n_p*t_c+t_r);
        
        while(pat_ptr<end_ptr){
            REAL_TYPE* end_col_ptr=pat_ptr+t_r-f_r+1;
            REAL_TYPE* start_next_col_ptr=pat_ptr+n_p;
            for(;pat_ptr<end_col_ptr;pat_ptr++){
                (*pat_ptr)=(*pat_ptr)+c;
            }
            pat_ptr=start_next_col_ptr;
        }
    }
    /*   printf("returning\n");*/
}


void direct_l2(const REAL_TYPE* img,int im_rows,int im_cols,const REAL_TYPE* pat,int n_p,int m_p,int* i_loc,int* j_loc,REAL_TYPE* vals,int ind_num){

    int ci,cj;
    const int* end_i=i_loc+ind_num;

    while(i_loc<end_i){
        cj=(*j_loc);
        
        while( (*j_loc)==cj && i_loc<end_i){
            ci=(*i_loc);
            (*vals)=l2_between_pat_and_win_in_img(img,pat,im_rows,im_cols,n_p,m_p,ci,cj);
            vals++;  
            i_loc++;
            j_loc++;
        }
    }
}


void direct_l2_2(const REAL_TYPE* img,int im_rows,int im_cols,const REAL_TYPE* pat,int n_p,int m_p,int* i_loc,int* j_loc,REAL_TYPE* vals,int ind_num){

    int cj;
    const int* end_i=i_loc+ind_num;
    const REAL_TYPE* patr =pat;
    const REAL_TYPE* winr;
    int* i_loc_ptr_to_firsr_in_row;
    int* j_loc_ptr_to_firsr_in_row;
    REAL_TYPE* vals_ptr_to_firsr_in_row;
    int pcolidx;
    REAL_TYPE diff;
    while(i_loc<end_i){
        cj=(*j_loc);

        i_loc_ptr_to_firsr_in_row=i_loc;
        j_loc_ptr_to_firsr_in_row=j_loc;
        vals_ptr_to_firsr_in_row=vals;
        *vals=0;


        for(pcolidx=0;pcolidx<m_p;pcolidx++){/*iterate the colums of pat*/

            if(pcolidx==0){
                *vals=0;
            }
            
            i_loc=i_loc_ptr_to_firsr_in_row;
            j_loc=j_loc_ptr_to_firsr_in_row;
            vals=vals_ptr_to_firsr_in_row;
        
            while( (*j_loc)==cj && i_loc<end_i){
                
                patr=pat+n_p*pcolidx;
                winr=(img+im_rows*cj+(*i_loc))+im_rows*pcolidx;
                
                while(patr<pat+n_p*(pcolidx+1)){

                    diff=(*patr)-(*winr);
                    (*vals) += (diff*diff);
                    patr++;
                    winr++;
                }
                i_loc++;
                j_loc++;
                vals++;
                
            }    
        }
    }
}

void direct_l2_3(const REAL_TYPE* img,int im_rows,int im_cols,const REAL_TYPE* pat,int n_p,int m_p,int* i_loc,int* j_loc,REAL_TYPE* vals,int ind_num){

    int cj;
    const int* end_i=i_loc+ind_num;
    const REAL_TYPE* patr =pat;
    const REAL_TYPE* winr;
    int* i_loc_ptr_to_firsr_in_row;
    int* j_loc_ptr_to_firsr_in_row;
    REAL_TYPE* vals_ptr_to_firsr_in_row;
    int pcolidx;
    REAL_TYPE diff;
    while(i_loc<end_i){
        cj=(*j_loc);

        i_loc_ptr_to_firsr_in_row=i_loc;
        j_loc_ptr_to_firsr_in_row=j_loc;
        vals_ptr_to_firsr_in_row=vals;
       
        while( (*j_loc)==cj && i_loc<end_i){

            *vals=0;
            patr=pat;
            winr=(img+im_rows*cj+(*i_loc));
            
            while(patr<pat+n_p){
                
                diff=(*patr)-(*winr);
                (*vals) += (diff*diff);
                patr++;
                winr++;
            }
            i_loc++;
            j_loc++;
            vals++;
            
        }
        
        for(pcolidx=1;pcolidx<m_p;pcolidx++){/*iterate the colums of pat*/

            
            i_loc=i_loc_ptr_to_firsr_in_row;
            j_loc=j_loc_ptr_to_firsr_in_row;
            vals=vals_ptr_to_firsr_in_row;
        
            while( (*j_loc)==cj && i_loc<end_i){
                
                patr=pat+n_p*pcolidx;
                winr=(img+im_rows*cj+(*i_loc))+im_rows*pcolidx;
                
                while(patr<pat+n_p*(pcolidx+1)){

                    diff=(*patr)-(*winr);
                    (*vals) += (diff*diff);
                    patr++;
                    winr++;
                }
                i_loc++;
                j_loc++;
                vals++;
                
            }    
        }
    }
}


REAL_TYPE inner_product(const REAL_TYPE* img,const REAL_TYPE* pat,int im_rows,int im_cols,int n_p,int m_p,int ci,int cj){
    
    REAL_TYPE p=0;
    REAL_TYPE m;
    const REAL_TYPE* pat_end_ptr=pat+(n_p*m_p);
    const REAL_TYPE* img_ptr=img+im_rows*cj+ci;
    const REAL_TYPE* img_end_col_ptr;
        
     while(pat<pat_end_ptr){
         img_end_col_ptr=img_ptr+n_p;
         for(;img_ptr<img_end_col_ptr;img_ptr++){
             m=(*pat)*(*img_ptr);
             p+=m;
             pat++;
    
         }
         img_ptr=img_ptr+im_rows-n_p;
     }
     return p;
}

/*
REAL_TYPE lidx_inner_product(const REAL_TYPE* img,const REAL_TYPE* pat,int im_rows,int im_cols,int n_p,int m_p,int lidx){
    
    REAL_TYPE p=0;
    REAL_TYPE m;
    const REAL_TYPE* pat_end_ptr=pat+(n_p*m_p);
    const REAL_TYPE* img_ptr=img+im_rows*cj+ci;
    const REAL_TYPE* img_end_col_ptr;
        
     while(pat<pat_end_ptr){
         img_end_col_ptr=img_ptr+n_p;
         for(;img_ptr<img_end_col_ptr;img_ptr++){
             m=(*pat)*(*img_ptr);
             p+=m;
             pat++;
    
         }
         img_ptr=img_ptr+im_rows-n_p;
     }
     return p;
}
*/
REAL_TYPE l2_between_pat_and_win_in_img(const REAL_TYPE* img,const REAL_TYPE* pat,int im_rows,int im_cols,int n_p,int m_p,int ci,int cj){
    
    REAL_TYPE p=0;
    REAL_TYPE diff;
    const REAL_TYPE* pat_end_ptr=pat+(n_p*m_p);
    const REAL_TYPE* img_ptr=img+im_rows*cj+ci;
    const REAL_TYPE* img_end_col_ptr;
        
     while(pat<pat_end_ptr){
         img_end_col_ptr=img_ptr+n_p;
         for(;img_ptr<img_end_col_ptr;img_ptr++){
             diff=(*pat)-(*img_ptr);
             p+=(diff*diff);
             pat++;
    
         }
         img_ptr=img_ptr+im_rows-n_p;
     }
     return p;
}



int U_before_loop_and_find(const REAL_TYPE* I,int I_r,REAL_TYPE* U,int from_r,int to_r,int from_c,int to_c, REAL_TYPE w,const REAL_TYPE* norms_of_wins_minus_mu,int img_r,int img_c,int n_p,int m_p,const REAL_TYPE*  wins_mu,int* i_loc,int* j_loc,REAL_TYPE threshold){

    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    
    const REAL_TYPE* I_br_ptr;
    const REAL_TYPE* I_bl_ptr;
    const REAL_TYPE* I_ur_ptr;
    const REAL_TYPE* I_ul_ptr;
     
    const REAL_TYPE* ul_const_term = I+from_r + (I_r*from_c);
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);

    /*
    REAL_TYPE* U_ptr=U;
    const REAL_TYPE* wins_mu_ptr = wins_mu;
    const REAL_TYPE* norms_of_wins_minus_mu_ptr =norms_of_wins_minus_mu;
    */
    
    REAL_TYPE tmp;

    int I_idx;
    int ci=0;
    int cj=0;
    REAL_TYPE* U_end = U+c_c*c_r;
    int support_sz=(to_r-from_r+1)*(to_c-from_c+1);
    
    int* i_next_ptr=i_loc;
    int* j_next_ptr=j_loc;
    
    while(U<U_end){
       
        I_idx = I_r*cj+ci;
      
        I_ul_ptr = ul_const_term+I_idx;
        I_ur_ptr = ur_const_term+I_idx;
        I_bl_ptr = bl_const_term+I_idx;
        I_br_ptr = br_const_term+I_idx;

        tmp = (*I_br_ptr);
        tmp -=(*I_ur_ptr);
        tmp -=(*I_bl_ptr);
        tmp +=(*I_ul_ptr);    
        tmp = w*(tmp-support_sz*(*wins_mu))/(*norms_of_wins_minus_mu);

        if(tmp>threshold){
            *U=tmp;
            (*i_next_ptr)=ci;
            (*j_next_ptr)=cj;
            i_next_ptr++;
            j_next_ptr++;
            /*num++;*/
        }
        /*
        else{
            *U_ptr=0;
        }
        */
        U++;
        wins_mu++;
        norms_of_wins_minus_mu++;
        ci++;
        if(ci==c_r){
            ci=0;
            cj++;
        }
    }
    /*return num;*/
    /*printf("%d \n",i_next_ptr-i_loc);*/
      return i_next_ptr-i_loc;
}  

int linear_idx_U_before_loop_and_find(const REAL_TYPE* I,int I_r,REAL_TYPE* U,int from_r,int to_r,int from_c,int to_c, REAL_TYPE w,const REAL_TYPE* norms_of_wins_minus_mu,int img_r,int img_c,int n_p,int m_p,const REAL_TYPE*  wins_mu,int* lidx,REAL_TYPE threshold){

    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    
    const REAL_TYPE* I_br_ptr;
    const REAL_TYPE* I_bl_ptr;
    const REAL_TYPE* I_ur_ptr;
    const REAL_TYPE* I_ul_ptr;
     
    const REAL_TYPE* ul_const_term = I+from_r + (I_r*from_c);
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);

    
    REAL_TYPE* U_ptr=U;
    /*
    const REAL_TYPE* wins_mu_ptr = wins_mu;
    const REAL_TYPE* norms_of_wins_minus_mu_ptr =norms_of_wins_minus_mu;
    */
    REAL_TYPE tmp;

    int I_idx;
    
    int ci=0;
    /*    int cj=0;*/
    
    REAL_TYPE* U_end = U+c_c*c_r;
    int support_sz=(to_r-from_r+1)*(to_c-from_c+1);
    
    int* lidx_ptr=lidx;
    /*int* lIidx_ptr=lIidx;*/

    int I_idx_col=0;
    while(U_ptr<U_end){

        /*        I_idx = I_r*cj+ci;*/
        I_idx = I_idx_col+ci;
      
        I_ul_ptr = ul_const_term+I_idx;
        I_ur_ptr = ur_const_term+I_idx;
        I_bl_ptr = bl_const_term+I_idx;
        I_br_ptr = br_const_term+I_idx;

        tmp = (*I_br_ptr);
        tmp -=(*I_ur_ptr);
        tmp -=(*I_bl_ptr);
        tmp +=(*I_ul_ptr);    
        tmp = w*(tmp-support_sz*(*wins_mu))/(*norms_of_wins_minus_mu);

        if(tmp>threshold){
            *U_ptr=tmp;
            (*lidx_ptr)=U_ptr-U;
            lidx_ptr++;
        }
        U_ptr++;
        wins_mu++;
        norms_of_wins_minus_mu++;
        ci++;

        /*
        if(ci==c_r){
            ci=0;
            cj++;
        }
        */
        if(ci==c_r){
            ci=0;
            I_idx_col+=I_r;
        }     
    }
    return lidx_ptr-lidx;
}  


/*
 indices_num = ncc_inside_loop(I,U,f_r,t_r,f_c,t_c,w_arr[t],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,i_loc2,j_loc2,i_loc,j_loc,threshold*residual_pat_norms[0]-residual_pat_norms[t]);
 }*/



int ncc_inside_loop(const REAL_TYPE* I,int I_r,REAL_TYPE* U, const REAL_TYPE* norms_of_wins_minus_mu,int c_r,int c_c,int from_r ,int to_r ,int from_c ,int to_c, REAL_TYPE w, const REAL_TYPE* wins_mu,int* i_ptr,int* j_ptr,int ind_num,REAL_TYPE threshold){

    int cj;

    /*
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    */
    const int* end_i_ptr=i_ptr+ind_num;
    int* i_ptr_for_next=i_ptr;
    int* j_ptr_for_next=j_ptr;
    
    const REAL_TYPE* I_br_ptr;
    const REAL_TYPE* I_bl_ptr;
    const REAL_TYPE* I_ur_ptr;
    const REAL_TYPE* I_ul_ptr;
    const REAL_TYPE* ul_const_term = I+from_r + I_r*from_c;
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);

    REAL_TYPE tmp;
    int I_col_term,I_idx;
    int C_col_term,C_idx;
    int support_sz=(to_r-from_r+1)*(to_c-from_c+1);

    while(i_ptr<end_i_ptr){

        cj=(*j_ptr);
        I_col_term = I_r*cj;
        C_col_term = c_r*cj;
        
        while( (*j_ptr)==cj && i_ptr<end_i_ptr){
            /*            ci=(*i_ptr);*/

            I_idx=I_col_term+(*i_ptr);
            C_idx=C_col_term+(*i_ptr);
            
            I_ul_ptr = ul_const_term+I_idx;
            I_ur_ptr = ur_const_term+I_idx;
            I_bl_ptr = bl_const_term+I_idx;
            I_br_ptr = br_const_term+I_idx;
            
            tmp = (*I_br_ptr);
            tmp -=(*I_ur_ptr);
            tmp -=(*I_bl_ptr);
            tmp +=(*I_ul_ptr);

            tmp=(w*tmp-w*support_sz*(*(wins_mu+C_idx)))/(*(norms_of_wins_minus_mu+C_idx));
            tmp=tmp+(*(U+C_idx));

            if(tmp>threshold){
                      
                (*(U+C_idx))=tmp;            
                (*i_ptr_for_next)=(*i_ptr);
                (*j_ptr_for_next)=cj;
                i_ptr_for_next++;
                j_ptr_for_next++;
                /*num++;*/
            }
            i_ptr++;
            j_ptr++;
        }
    }
    /*printf("%d \n",num);*/
    /*    return num;*/
    return ind_num-(end_i_ptr-i_ptr_for_next);
}


void direct_conv(const REAL_TYPE* img,int im_rows,int im_cols,const REAL_TYPE* pat,int n_p,int m_p,int* i_loc,int* j_loc,REAL_TYPE* vals,int ind_num){

    int ci,cj;
    const int* end_i=i_loc+ind_num;

    while(i_loc<end_i){
        cj=(*j_loc);
        
        while( (*j_loc)==cj && i_loc<end_i){
            ci=(*i_loc);
            (*vals)=inner_product(img,pat,im_rows,im_cols,n_p,m_p,ci,cj);
            vals++;  
            i_loc++;
            j_loc++;
        }
    }
}

int lidx_ncc_inside_loop(const REAL_TYPE* I,int I_r,REAL_TYPE* U, const REAL_TYPE* norms_of_wins_minus_mu,int c_r,int c_c,int from_r ,int to_r ,int from_c ,int to_c, REAL_TYPE w, const REAL_TYPE* wins_mu,int* lidx,int ind_num,REAL_TYPE threshold){

    
    int ci,cj;

    const int* end_lidx=lidx+ind_num;

    int* lidx_ptr_for_next=lidx;
    
    const REAL_TYPE* I_br_ptr;
    const REAL_TYPE* I_bl_ptr;
    const REAL_TYPE* I_ur_ptr;
    const REAL_TYPE* I_ul_ptr;
    const REAL_TYPE* ul_const_term = I+from_r + I_r*from_c;
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);

    REAL_TYPE tmp;
    int I_idx;
    int support_sz=(to_r-from_r+1)*(to_c-from_c+1);

     cj=(*lidx)/c_r;
     ci=(*lidx)-c_r*cj;

     int I_idx_col=I_r*cj;
    while(lidx<end_lidx){

       
        I_idx = I_idx_col+ci;        
        I_ul_ptr = ul_const_term+I_idx;
        I_ur_ptr = ur_const_term+I_idx;
        I_bl_ptr = bl_const_term+I_idx;
        I_br_ptr = br_const_term+I_idx;
        
        tmp = (*I_br_ptr);
        tmp -=(*I_ur_ptr);
        tmp -=(*I_bl_ptr);
        tmp +=(*I_ul_ptr);
        
        tmp=(w*tmp-w*support_sz*(*(wins_mu+(*lidx))))/(*(norms_of_wins_minus_mu+(*lidx)));
        tmp=tmp+(*(U+(*lidx)));
        
        if(tmp>threshold){
            
            (*(U+(*lidx)))=tmp;            
            (*lidx_ptr_for_next)=(*lidx);
            lidx_ptr_for_next++;
            /*num++;*/
        }
        lidx++;

        ci=ci+(*lidx)-(*(lidx-1));
        if(ci>=c_r){
            cj=(*lidx)/c_r;
            ci=(*lidx)-c_r*cj;
            I_idx_col=I_r*cj;
        }
            
    }
    /*printf("%d \n",num);*/
    return ind_num-(end_lidx-lidx_ptr_for_next);
}


void lidx_direct_ncc(const REAL_TYPE* img,int im_rows,int im_cols,const REAL_TYPE* pat,int n_p,int m_p,int* lidx,REAL_TYPE* vals,int ind_num,REAL_TYPE pat_norm,const REAL_TYPE* norms_of_wins_minus_mu){

    int ci,cj;
    const int* end_lidx=lidx+ind_num;
    REAL_TYPE d;
    int c_r=im_rows-n_p+1;


    while(lidx<end_lidx){

        ci=(*lidx)%c_r;
        cj=((*lidx)-ci)/c_r;
    
        d=pat_norm*(*(norms_of_wins_minus_mu+ci+c_r*cj));
        if(d>0.001){
            (*vals)=inner_product(img,pat,im_rows,im_cols,n_p,m_p,ci,cj)/d;
        }
        else{
            *vals=0;
        }

        vals++;  
        lidx++;
    }
}



void direct_ncc(const REAL_TYPE* img,int im_rows,int im_cols,const REAL_TYPE* pat,int n_p,int m_p,int* i_loc,int* j_loc,REAL_TYPE* vals,int ind_num,REAL_TYPE pat_norm,const REAL_TYPE* norms_of_wins_minus_mu){

    int ci,cj;
    const int* end_i=i_loc+ind_num;
    REAL_TYPE d;
    int c_r=im_rows-n_p+1;
    while(i_loc<end_i){
        cj=(*j_loc);
        
        while( (*j_loc)==cj && i_loc<end_i){
            ci=(*i_loc);
            d=pat_norm*(*(norms_of_wins_minus_mu+ci+c_r*cj));
            if(d>0.001){
                (*vals)=inner_product(img,pat,im_rows,im_cols,n_p,m_p,ci,cj)/d;
            }
            else{
                *vals=0;
            }
            /*
            if(ci>94 && ci<97 && (cj>466 & cj<470)){
                printf(" d=%f cor=%f\n",d,*vals);
            }
            
            (*vals)=(*vals)/d;
            */
            vals++;  
            i_loc++;
            j_loc++;
        }
    }
}


void wins_mu_and_norms_of_wins_minus_mu(const REAL_TYPE* I,const REAL_TYPE* I2, int img_r,int img_c, int n_p,int m_p, REAL_TYPE* wins_mu, REAL_TYPE* norms_of_wins_minus_mu){    


    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int I_r=img_r+1;

    /*printf("c_r=%d, c_c=%d\n",c_r,c_c);*/

    int win_size=n_p*m_p;
    
    REAL_TYPE* C_ptr=wins_mu;
    REAL_TYPE* C_end_ptr=wins_mu+c_r*c_c;
    REAL_TYPE* norms_ptr=norms_of_wins_minus_mu;
    
    const REAL_TYPE* I_br_ptr=I+I_r*m_p+n_p;
    const REAL_TYPE* I_bl_ptr=I+n_p;
    const REAL_TYPE* I_ur_ptr=I+I_r*m_p;
    const REAL_TYPE* I_ul_ptr=I;

    const REAL_TYPE* I2_br_ptr=I2+I_r*m_p+n_p;
    const REAL_TYPE* I2_bl_ptr=I2+n_p;
    const REAL_TYPE* I2_ur_ptr=I2+I_r*m_p;
    const REAL_TYPE* I2_ul_ptr=I2;

    while(C_ptr<C_end_ptr){
        REAL_TYPE* next_col_ptr = C_ptr+c_r;
        while(C_ptr<next_col_ptr){

            REAL_TYPE tmp1 = (*I_br_ptr);
            tmp1 -=(*I_ur_ptr);
            tmp1 -=(*I_bl_ptr);
            tmp1 +=(*I_ul_ptr);
            (*C_ptr)=tmp1/win_size;


            REAL_TYPE tmp2 = (*I2_br_ptr);
            tmp2 -=(*I2_ur_ptr);
            tmp2 -=(*I2_bl_ptr);
            tmp2 +=(*I2_ul_ptr);
            (*norms_ptr)=sqrt(tmp2-((*C_ptr)*tmp1));

            
            I_br_ptr++;
            I_ur_ptr++;
            I_bl_ptr++;
            I_ul_ptr++;
            C_ptr++;
            I2_br_ptr++;
            I2_ur_ptr++;
            I2_bl_ptr++;
            I2_ul_ptr++;
            norms_ptr++;
        }

        I_br_ptr=I_br_ptr+n_p;
        I_ur_ptr=I_ur_ptr+n_p;
        I_bl_ptr=I_bl_ptr+n_p;
        I_ul_ptr=I_ul_ptr+n_p;

        I2_br_ptr=I2_br_ptr+n_p;
        I2_ur_ptr=I2_ur_ptr+n_p;
        I2_bl_ptr=I2_bl_ptr+n_p;
        I2_ul_ptr=I2_ul_ptr+n_p;

    }
}

void wins_mu_and_norms_of_wins_minus_mu2(const REAL_TYPE* I,const REAL_TYPE* I2, int img_r,int img_c, int n_p,int m_p, REAL_TYPE* wins_mu, REAL_TYPE* norms_of_wins_minus_mu){    

   
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int I_r=img_r+1;

    int win_size=n_p*m_p;
    
    REAL_TYPE* C_ptr=wins_mu;
    REAL_TYPE* C_end_ptr=wins_mu+c_r*c_c;
    REAL_TYPE* norms_ptr=norms_of_wins_minus_mu;
    
    const REAL_TYPE* I_br_ptr0;
    const REAL_TYPE* I_bl_ptr0;
    const REAL_TYPE* I_ur_ptr0;
    const REAL_TYPE* I_ul_ptr0;
    
    const REAL_TYPE* I2_br_ptr;
    const REAL_TYPE* I2_bl_ptr;
    const REAL_TYPE* I2_ur_ptr;
    const REAL_TYPE* I2_ul_ptr;

    int ur_const_term0 = I_r*m_p;
    int br_const_term0 = n_p + I_r*m_p;

    int I_idx_col =0;
    int I_idx,ci=0;
    while(C_ptr<C_end_ptr){
        I_idx = I_idx_col+ci;

        I_ul_ptr0 = I+I_idx;
        I_ur_ptr0 = I+ur_const_term0+I_idx;
        I_bl_ptr0 = I+n_p+I_idx;
        I_br_ptr0 = I+br_const_term0+I_idx;

        I2_ul_ptr = I2+I_idx;
        I2_ur_ptr = I2+ur_const_term0+I_idx;
        I2_bl_ptr = I2+n_p+I_idx;
        I2_br_ptr = I2+br_const_term0+I_idx;

        REAL_TYPE tmp1 = (*I_br_ptr0);
        tmp1 -=(*I_ur_ptr0);
        tmp1 -=(*I_bl_ptr0);
        tmp1 +=(*I_ul_ptr0);
        (*C_ptr)=tmp1/win_size;
        REAL_TYPE tmp2 = (*I2_br_ptr);
        tmp2 -=(*I2_ur_ptr);
        tmp2 -=(*I2_bl_ptr);
        tmp2 +=(*I2_ul_ptr);
        (*norms_ptr)=sqrt(tmp2-((*C_ptr)*tmp1));
        C_ptr++;
        norms_ptr++;
        ci++;
        if(ci==c_r){
            ci=0;
            I_idx_col+=I_r;
        }  
        
    }
    
}



/*
wins_sums= box_corr(I,[1,n_p,1,m_p],1,n_p,m_p);
wins_mu=wins_sums/(n_p*m_p);
wins_norm2 = box_corr(I2,[1,n_p,1,m_p],1,n_p,m_p);
norms_of_wins_minus_mu = sqrt(abs(wins_norm2-(wins_mu.*wins_mu*n_p*m_p)));
*/

/*int indexes_num = linear_idx_norms_mu_U_before_loop_find(I,img_r+1,U,box_arr[0]-1,box_arr[1]-1,box_arr[2]-1,box_arr[3]-1,w_arr[0],norms_of_wins_minus_mu,img_r,img_c,n_p,m_p, wins_mu,lidx,threshold*residual_pat_norms[0]-residual_pat_norms[1]);
 */
int linear_idx_norms_mu_U_before_loop_find(const REAL_TYPE* I,const REAL_TYPE* I2,int I_r,REAL_TYPE* U,int from_r,int to_r,int from_c,int to_c, REAL_TYPE w,REAL_TYPE* norms_of_wins_minus_mu,int img_r,int img_c,int n_p,int m_p,REAL_TYPE*  wins_mu,int* lidx,REAL_TYPE threshold){

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int win_size=n_p*m_p;
    int ur_const_term0 = I_r*m_p;
    int br_const_term0 = n_p + I_r*m_p;

    const REAL_TYPE* I_ul_ptr0=I;
    const REAL_TYPE* I_bl_ptr0=I_ul_ptr0+n_p;
    const REAL_TYPE* I_ur_ptr0=I_ul_ptr0+ur_const_term0;
    const REAL_TYPE* I_br_ptr0=I_ul_ptr0+br_const_term0;    

    const REAL_TYPE* I2_ul_ptr=I2;
    const REAL_TYPE* I2_bl_ptr=I2_ul_ptr+n_p;
    const REAL_TYPE* I2_ur_ptr=I2_ul_ptr+ur_const_term0;
    const REAL_TYPE* I2_br_ptr=I2_ul_ptr+br_const_term0;    

    const REAL_TYPE* ul_const_term = I+from_r + (I_r*from_c);
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);

    const REAL_TYPE* I_ul_ptr = ul_const_term;
    const REAL_TYPE* I_ur_ptr = ur_const_term;
    const REAL_TYPE* I_bl_ptr = bl_const_term;
    const REAL_TYPE* I_br_ptr = br_const_term;

    REAL_TYPE* U_ptr=U;
   
    REAL_TYPE tmp,tmp_wins_mu,tmp_norms;
    int ci=0;    
    REAL_TYPE* U_end = U+c_c*c_r;
    int support_sz=(to_r-from_r+1)*(to_c-from_c+1);
    int* lidx_ptr=lidx;
    int n_p_plus_one=n_p+1;

    while(U_ptr<U_end){
        tmp_wins_mu = (*I_ul_ptr0);
        tmp_wins_mu -=(*I_ur_ptr0);
        tmp_wins_mu -=(*I_bl_ptr0);
        tmp_wins_mu +=(*I_br_ptr0);
        (*wins_mu)=tmp_wins_mu/win_size;

        tmp_norms = (*I2_br_ptr);
        tmp_norms -=(*I2_ur_ptr);
        tmp_norms -=(*I2_bl_ptr);
        tmp_norms +=(*I2_ul_ptr);
        (*norms_of_wins_minus_mu)=sqrt(tmp_norms-((*wins_mu)*tmp_wins_mu));
              
        tmp = (*I_br_ptr);
        tmp -=(*I_ur_ptr);
        tmp -=(*I_bl_ptr);
        tmp +=(*I_ul_ptr);    
        tmp = w*(tmp-support_sz*(*wins_mu))/(*norms_of_wins_minus_mu);

        /*** TODO - REMOVE THIS LINE! ***/
        /**U_ptr=tmp;*/
        if(tmp>threshold){
            *U_ptr=tmp;
            (*lidx_ptr)=U_ptr-U;
            lidx_ptr++;
        }
        U_ptr++;
        wins_mu++;
        norms_of_wins_minus_mu++;
        ci++;
        if(ci<c_r){
            I_ul_ptr0++;
            I_ur_ptr0++;
            I_bl_ptr0++;
            I_br_ptr0++;
            I2_ul_ptr++;
            I2_ur_ptr++;
            I2_bl_ptr++;
            I2_br_ptr++;
            I_ul_ptr++;
            I_ur_ptr++;
            I_bl_ptr++;
            I_br_ptr++;
        }
        else{
            ci=0;
            I_ul_ptr0+=n_p_plus_one;
            I_ur_ptr0+=n_p_plus_one;
            I_bl_ptr0+=n_p_plus_one;
            I_br_ptr0+=n_p_plus_one;

            I2_ul_ptr+=n_p_plus_one;
            I2_ur_ptr+=n_p_plus_one;
            I2_bl_ptr+=n_p_plus_one;
            I2_br_ptr+=n_p_plus_one;
            
            I_ul_ptr+=n_p_plus_one;
            I_ur_ptr+=n_p_plus_one;
            I_bl_ptr+=n_p_plus_one;
            I_br_ptr+=n_p_plus_one;
        }  
    }
    return lidx_ptr-lidx;
}


/**********************************************************************************/

int lidx_ncc_inside_loop_scan_all_img(const REAL_TYPE* I,int I_r,REAL_TYPE* U, const REAL_TYPE* norms_of_wins_minus_mu,int c_r,int c_c,int from_r ,int to_r ,int from_c ,int to_c, REAL_TYPE w, const REAL_TYPE* wins_mu,int* lidx,int ind_num,REAL_TYPE threshold,int n_p){

    const REAL_TYPE* ul_const_term = I+from_r + (I_r*from_c);
    const REAL_TYPE* ur_const_term = I+from_r + I_r*(1+to_c);
    const REAL_TYPE* bl_const_term = I+1+to_r + I_r*from_c;
    const REAL_TYPE* br_const_term = I+1+to_r + I_r*(1+to_c);
    const REAL_TYPE* I_ul_ptr = ul_const_term;
    const REAL_TYPE* I_ur_ptr = ur_const_term;
    const REAL_TYPE* I_bl_ptr = bl_const_term;
    const REAL_TYPE* I_br_ptr = br_const_term;

    REAL_TYPE tmp;
    int ci=0;    
    int support_sz=(to_r-from_r+1)*(to_c-from_c+1);
    int n_p_plus_one=n_p+1;
    int current=0;
    int* lidx_ptr_for_next=lidx;
    const int* end_lidx=lidx+ind_num;
    while(lidx<end_lidx){
        if(current==(*lidx)){
            tmp = (*I_br_ptr);
            tmp -=(*I_ur_ptr);
            tmp -=(*I_bl_ptr);
            tmp +=(*I_ul_ptr);    
            tmp = w*(tmp-support_sz*(*wins_mu))/(*norms_of_wins_minus_mu);
            tmp=tmp+(*U);
            
            if(tmp>threshold){
                *U=tmp;
                (*lidx_ptr_for_next)=(*lidx);
                lidx_ptr_for_next++;
            }
            lidx++;
        }
        current++;
        U++;
        wins_mu++;
        norms_of_wins_minus_mu++;
        ci++;
        if(ci<c_r){
            I_ul_ptr++;
            I_ur_ptr++;
            I_bl_ptr++;
            I_br_ptr++;
        }
        else{
            ci=0;
            I_ul_ptr+=n_p_plus_one;
            I_ur_ptr+=n_p_plus_one;
            I_bl_ptr+=n_p_plus_one;
            I_br_ptr+=n_p_plus_one;
        }  
    }   
    return ind_num-(end_lidx-lidx_ptr_for_next);
}
/***********************************************************************************/
void linear_idx_norms_mu(const REAL_TYPE* I,const REAL_TYPE* I2,int I_r,REAL_TYPE* norms_of_wins_minus_mu,int img_r,int img_c,int n_p,int m_p,REAL_TYPE*  wins_mu){

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    int win_size=n_p*m_p;
    int ur_const_term0 = I_r*m_p;
    int br_const_term0 = n_p + I_r*m_p;

    const REAL_TYPE* I_ul_ptr0=I;
    const REAL_TYPE* I_bl_ptr0=I_ul_ptr0+n_p;
    const REAL_TYPE* I_ur_ptr0=I_ul_ptr0+ur_const_term0;
    const REAL_TYPE* I_br_ptr0=I_ul_ptr0+br_const_term0;    

    const REAL_TYPE* I2_ul_ptr=I2;
    const REAL_TYPE* I2_bl_ptr=I2_ul_ptr+n_p;
    const REAL_TYPE* I2_ur_ptr=I2_ul_ptr+ur_const_term0;
    const REAL_TYPE* I2_br_ptr=I2_ul_ptr+br_const_term0;    

    REAL_TYPE tmp_wins_mu,tmp_norms;
    int ci=0;        
    int n_p_plus_one=n_p+1;
    REAL_TYPE* wins_mu_end=wins_mu+c_r*c_c;
    while(wins_mu<wins_mu_end){
        tmp_wins_mu = (*I_ul_ptr0);
        tmp_wins_mu -=(*I_ur_ptr0);
        tmp_wins_mu -=(*I_bl_ptr0);
        tmp_wins_mu +=(*I_br_ptr0);
        (*wins_mu)=tmp_wins_mu/win_size;

        tmp_norms = (*I2_br_ptr);
        tmp_norms -=(*I2_ur_ptr);
        tmp_norms -=(*I2_bl_ptr);
        tmp_norms +=(*I2_ul_ptr);
        (*norms_of_wins_minus_mu)=sqrt(tmp_norms-((*wins_mu)*tmp_wins_mu));
              
        wins_mu++;
        norms_of_wins_minus_mu++;
        ci++;
        if(ci<c_r){
            I_ul_ptr0++;
            I_ur_ptr0++;
            I_bl_ptr0++;
            I_br_ptr0++;
            I2_ul_ptr++;
            I2_ur_ptr++;
            I2_bl_ptr++;
            I2_br_ptr++;
        }
        else{
            ci=0;
            I_ul_ptr0+=n_p_plus_one;
            I_ur_ptr0+=n_p_plus_one;
            I_bl_ptr0+=n_p_plus_one;
            I_br_ptr0+=n_p_plus_one;

            I2_ul_ptr+=n_p_plus_one;
            I2_ur_ptr+=n_p_plus_one;
            I2_bl_ptr+=n_p_plus_one;
            I2_br_ptr+=n_p_plus_one;
        }  
    }
}


/**********************************************************************************/
