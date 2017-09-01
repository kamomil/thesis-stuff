#include <stdlib.h>
/*#include <sys/types.h>*/
#include <stdio.h>
/*#include <string.h>*/
#include <sys/time.h>



#define SEC_RES 1
#define IDX_LENGTH 100000


void mod_vs_dev(){

    int i;
    int n[3]={1,2,3};
    int idx[IDX_LENGTH];
    int s;
    idx[0]=10;
    for(i=1;i<IDX_LENGTH;i++){
        idx[i]=idx[i-1]+n[i%3];
        //        printf("%d ",idx[i]);
    }
    
    struct timeval  first_tv1, first_tv2;
    struct timeval  sec_tv1, sec_tv2;

    
    gettimeofday(&sec_tv1,NULL);
    s=idx[0];
    for(i=1;i<IDX_LENGTH;i++){
        s=idx[i]/100;
        printf("%d ",s);
    }
    gettimeofday(&sec_tv2,NULL);

    
    gettimeofday(&first_tv1,NULL);
    for(i=0;i<IDX_LENGTH;i++){
        s=idx[i]%100;
    }
    i=s;
    gettimeofday(&first_tv2,NULL);

    double first_run_time=SEC_RES * ((double) (first_tv2.tv_usec - first_tv1.tv_usec) / 1000000 +(double) (first_tv2.tv_sec - first_tv1.tv_sec));

    double sec_run_time=SEC_RES * ((double) (sec_tv2.tv_usec - sec_tv1.tv_usec) / 1000000 +(double) (sec_tv2.tv_sec - sec_tv1.tv_sec));

    printf("mod %f, dev %f\n",first_run_time, sec_run_time);



}



int main(){

    mod_vs_dev();
    int i;
    int n[3]={1,1,1};
    int idx[IDX_LENGTH];

    idx[0]=10;
    for(i=1;i<IDX_LENGTH;i++){
        idx[i]=idx[i-1]+n[i%3];
        //        printf("%d ",idx[i]);
    }
    
    int s;

    struct timeval  first_tv1, first_tv2;
    struct timeval  sec_tv1, sec_tv2;

    
    gettimeofday(&sec_tv1,NULL);
    s=idx[0];
    for(i=1;i<IDX_LENGTH;i++){
        s++;//(idx[i]-idx[i-1]);
    }
    gettimeofday(&sec_tv2,NULL);

    
    gettimeofday(&first_tv1,NULL);
    for(i=0;i<IDX_LENGTH;i++){
        s=idx[i];
    }
    gettimeofday(&first_tv2,NULL);

    double first_run_time=SEC_RES * ((double) (first_tv2.tv_usec - first_tv1.tv_usec) / 1000000 +(double) (first_tv2.tv_sec - first_tv1.tv_sec));

    double sec_run_time=SEC_RES * ((double) (sec_tv2.tv_usec - sec_tv1.tv_usec) / 1000000 +(double) (sec_tv2.tv_sec - sec_tv1.tv_sec));

    printf("first %f, sec %f\n",first_run_time, sec_run_time);
    return 0;
}
