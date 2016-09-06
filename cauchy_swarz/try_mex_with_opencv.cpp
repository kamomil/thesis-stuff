#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>

using namespace cv;

#define MAX_IMG_SIZE 1000000
#define MAX_BOX_NUM 3000

/* prhs - the array of inputs (each cell contain the pointer and size), will */
/* also have pointer to the output */
void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* Pointer to an mxArray of type float */
    const float* I = (const float*)mxGetPr(prhs[0]); 
    int img_r = mxGetM(prhs[0]); 
    int img_c = mxGetN(prhs[0]); 

    float* I2 = (float*) mxMalloc(sizeof(float)*img_r*img_c);

    int i;
    for(i=0;i<img_r*img_c;i++){
        I2[i]=I[i];
    }

    const float* reconst_pat = (const float*)mxGetPr(prhs[1]); 
    int n_p = mxGetM(prhs[4]); 
    int m_p = mxGetN(prhs[4]); 

    float* pat2 =  (float*)  mxMalloc(sizeof(float)*n_p*m_p);

    for(i=0;i<n_p*m_p;i++){
        pat2[i]=reconst_pat[i];
    }

    
    int c_r=img_r-n_p+1;
    int c_c=img_c-m_p+1;

    plhs[0]= mxCreateNumericMatrix(c_r,c_c,mxSINGLE_CLASS ,mxREAL);

    float* result_arr =  (float*)mxGetPr(plhs[0]);
    
    /*    const Mat img(5, 6, CV_32F,const_cast<float *>I);*/
    printf("1\n");
    Mat img(img_r, img_c, CV_32F,I2);
    printf("2\n");
    Mat templ(n_p, m_p, CV_32F,pat2);

    int result_cols = img.cols - templ.cols + 1;
    int result_rows = img.rows - templ.rows + 1;
    printf("3\n");
    Mat result( result_cols, result_rows, CV_32F );
    printf("4\n");
    matchTemplate(img, templ ,  result, CV_TM_SQDIFF);
    printf("5\n");
    /*  float* matData = (float*)myMat.data;*/
    result_arr = (float*)result.data;                   
}

