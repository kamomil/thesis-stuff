#include "mex.h"
#include "matrix.h"
#include <stdio.h>

#define MAX_IMG_SIZE 1000000
#define MAX_BOX_NUM 3000

/* REMEMBER a RxC matrix is represented such as [1, ... R, R+1, ....  2R ....] */
/* That is, C consecutive blocks each of size R, the ith block represent the ith column */

void box_corr_given_indices(const double* I, int img_r,int img_col, const int* box_arr,int box_num, const double* w_arr,int n_p,int m_p,const int* tpi,const int* tpj , int ind_num, double* C);

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


   
    if(mxGetClassID(prhs[1]) == mxDOUBLE_CLASS)
        mexErrMsgTxt("box_arr should be int (recieved double)");
    
    const int* box_arr = (int *) mxGetData(prhs[1]);
    int box_num = mxGetM(prhs[1]);

 
    const double* w_arr =  mxGetPr(prhs[2]);
    int n_p = (int)mxGetScalar(prhs[3]);
    int m_p = (int)mxGetScalar(prhs[4]); 

    if(mxGetClassID(prhs[5]) == mxDOUBLE_CLASS)
        mexErrMsgTxt("indicese should be int (recieved double)");
    
    if(mxGetClassID(prhs[6]) == mxDOUBLE_CLASS)
        mexErrMsgTxt("indices should be int (recieved double)");

    int* tpi = (int *) mxGetData(prhs[5]);
    int* tpj = (int *) mxGetData(prhs[6]);
    int ind_num = mxGetM(prhs[5]);
    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;
    
    
    plhs[0]=mxCreateDoubleMatrix(c_r,c_c,mxREAL);
    box_corr_given_indices(I,img_r,img_c,box_arr,box_num,w_arr,n_p,m_p,tpi,tpj,ind_num, mxGetPr(plhs[0]));
}

void box_corr_given_indices(const double* I,int img_r,int img_c, const int* box_arr,int box_num, const double* w_arr,int n_p,int m_p,const int* tpi,const int* tpj , int ind_num, double* C)
{
   
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    int bidx,idxidx;
    int ci,cj,Ij,Ii;

    /*    for(bidx=0;bidx<box_num;bidx++){
        int from_r = box_arr[bidx]-1;
        int to_r = box_arr[bidx+box_num]-1;
        int from_c = box_arr[bidx+2*box_num]-1;
        int to_c = box_arr[bidx+3*box_num]-1;
           
        for(idxidx=0;idxidx<ind_num;idxidx++){

            cj=tpi[idxidx]-1;
            ci=tpj[idxidx]-1; //columns , sorted by find, for fast column scan 

            Ij=cj+1;
            Ii=ci+1;
            
            double tmp = I[(Ij+to_r)+(img_r+1)*(Ii+to_c)];
            tmp=tmp-I[(Ij+to_r)+(img_r+1)*(Ii+from_c-1)];
            tmp=tmp-I[(Ij+from_r-1)+(img_r+1)*(Ii+to_c)];
            tmp=tmp+I[(Ij+from_r-1)+(img_r+1)*(Ii+from_c-1)];
            
            C[c_r*ci+cj]=C[c_r*ci+cj]+w_arr[bidx]*tmp;
        }     
    }*/

    /*printf("bo\n");*/
    for(bidx=0;bidx<box_num;bidx++){
        int from_r = box_arr[bidx]-1;
        int to_r = box_arr[bidx+box_num]-1;
        int from_c = box_arr[bidx+2*box_num]-1;
        int to_c = box_arr[bidx+3*box_num]-1;

        idxidx=0;
        cj=0;

        while(cj<c_c && idxidx<ind_num){
            cj=tpj[idxidx]-1;
            Ij=cj+1;            
            while(tpj[idxidx]-1==cj && idxidx<ind_num){

                ci=tpi[idxidx]-1;
                Ii=ci+1;
                    
                double tmp = I[(Ii+to_r)+(img_r+1)*(Ij+to_c)];
                tmp=tmp-I[(Ii+to_r)+(img_r+1)*(Ij+from_c-1)];
                tmp=tmp-I[(Ii+from_r-1)+(img_r+1)*(Ij+to_c)];
                tmp=tmp+I[(Ii+from_r-1)+(img_r+1)*(Ij+from_c-1)];
                
                C[c_r*cj+ci]=C[c_r*cj+ci]+w_arr[bidx]*tmp;
                idxidx++;
            }
            
        }     
    }

}
