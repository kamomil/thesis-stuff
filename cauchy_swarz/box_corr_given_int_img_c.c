#include "mex.h"
#include "matrix.h"
#include <stdio.h>

#define MAX_IMG_SIZE 1000000
#define MAX_BOX_NUM 3000

/* REMEMBER a RxC matrix is represented such as [1, ... R, R+1, ....  2R ....] */
/* That is, C consecutive blocks each of size R, the ith block represent the ith column */

void box_corr_given_int_img(const double* I, int img_r,int img_col, const int* box_arr,int box_num, const double* w_arr,int n_p,int m_p, double* C);

/* prhs - the array of inputs (each cell contain the pointer and size), will */
/* also have pointer to the output */
void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* Pointer to an mxArray of type double */
    const double* I = mxGetPr(prhs[0]); 
    int img_r = mxGetM(prhs[0])-1; 
    int img_c = mxGetN(prhs[0])-1; 


    mxClassID category = mxGetClassID(prhs[1]);

    if(category == mxDOUBLE_CLASS)
        mexErrMsgTxt("box_arr should be int (recieved double)");
    
    const int* box_arr = (int *) mxGetData(prhs[1]);
    int box_num = mxGetM(prhs[1]);

 
    const double* w_arr =  mxGetPr(prhs[2]);
    int n_p = (int)mxGetScalar(prhs[3]);
    int m_p = (int)mxGetScalar(prhs[4]); 

    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    
    plhs[0]=mxCreateDoubleMatrix(c_r,c_c,mxREAL);
    box_corr_given_int_img(I,img_r,img_c,box_arr,box_num,w_arr,n_p,m_p, mxGetPr(plhs[0]));
}

/*
void print_matrix(double* m,int r,int c){

    int i;
    int j;
    int k;
    for(i=0;i<r;i++){
        j=0;
        for(k=0;k<c;k++){
            printf("%f ",m[i+j]);
            j+=r;
        }
        printf("\n");
    }
}
*/

void box_corr_given_int_img(const double* I,int img_r,int img_c, const int* box_arr,int box_num, const double* w_arr,int n_p,int m_p, double* C)
{
   
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    int bidx=0;
    int ci,cj,Ii,Ij;

    for(bidx=0;bidx<box_num;bidx++){
        int from_r = box_arr[bidx]-1;
        int to_r = box_arr[bidx+box_num]-1;
        int from_c = box_arr[bidx+2*box_num]-1;
        int to_c = box_arr[bidx+3*box_num]-1;
   
        Ii=1;
        for(ci=0;ci<c_c;ci++){
            Ij=1;
            for(cj=0;cj<c_r;cj++){
            
                double tmp = I[(Ij+to_r)+(img_r+1)*(Ii+to_c)];
                tmp=tmp-I[(Ij+to_r)+(img_r+1)*(Ii+from_c-1)];
                tmp=tmp-I[(Ij+from_r-1)+(img_r+1)*(Ii+to_c)];
                tmp=tmp+I[(Ij+from_r-1)+(img_r+1)*(Ii+from_c-1)];
                
                C[c_r*ci+cj]=C[c_r*ci+cj]+w_arr[bidx]*tmp;

                Ij++;
            }
            Ii++;
        }     
    }
}
