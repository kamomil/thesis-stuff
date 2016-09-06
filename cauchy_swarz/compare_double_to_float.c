#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

void db(int n){
    double x[1000000];
    
    volatile int i=0;
    for (i=0; i<100000; i++){
        x[i]=((double ) n+i)/2;
    }

    for (i=3; i<100000; i++){
        x[i]=x[i-1]*x[i-3];
    }
}


void f(int n){
    float x[1000000];
    
    volatile int i=0;
    for (i=0; i<100000; i++){
        x[i]=((float ) n+i)/2;
    }

    for (i=3; i<100000; i++){
        x[i]=x[i-1]*x[i-3];
    }
}



int main(){

    printf("maximu integer is %d\n",INT_MAX);
    printf("maximu char is %d\n",CHAR_MAX);
    printf("maximu short is %d\n",SHRT_MAX);
    printf("size of float is %lu\n",sizeof(float));
    printf("size of double is %lu\n",sizeof(double));
    
    time_t start1,end1;
     double dif1;
     time (&start1);
     
    volatile int i=0;
    for (i=0; i<30000; i++){
        db(i);
    }

    time(&end1);
    dif1 = difftime (end1,start1);

 
    /////////////////////////////////////////////////////

     time_t start2,end2;
     double dif2;

     time(&start2);
     
    for (i=0; i<30000; i++){
        f(i);
    }

     time(&end2);
    dif2 = difftime (end2,start2);

  
  
  

    printf("time to db: %f\n",dif1);
    printf("time to f:  %f\n",dif2);

    return 0;
}
