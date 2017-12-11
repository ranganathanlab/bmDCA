#include <stdlib.h>
#include <stdio.h>
#include "math.h"

int main(int argc, char *argv[]){

int i,j,k,a,b,m;
int M,N,q;
double Meff;
FILE *fp;
fp=fopen(argv[1],"r");
//fp=fopen("../Model_inference/BetaLactamaseN10Q21/align.txt","r");
fscanf(fp,"%d %d %d\n",&M,&N,&q);

FILE *fp2;
fp2=fopen(argv[2],"r");
//fp2=fopen("../Model_inference/BetaLactamaseN10Q21/weights.txt","r");

int *x;
x=(int *)malloc(sizeof(int)*N*M);

double *n1;
n1=(double *)malloc(sizeof(double)*q);

double *n2;
n2=(double *)malloc(sizeof(double)*q*q);

double *n3;
n3=(double *)malloc(sizeof(double)*q*q*q);


double *W;
W=(double *)malloc(sizeof(double)*M);


//reading...
for(m=0;m<M;m++){
for(i=0;i<N;i++){
fscanf(fp,"%d",&x[m*N+i]);
}
fscanf(fp,"\n");
}

for(m=0;m<M;m++){
fscanf(fp2,"%lf\n",&W[m]);
}

//counting...

Meff=0;

for(m=0;m<M;m++){

//printf("%lf ",W[m]);
Meff+=W[m];
}


//output
FILE *fpw;
fpw=fopen("stat_align_1p.txt","w");
for(i=0;i<N;i++){
fprintf(fpw,"%d ",i);
for(a=0;a<q;a++){n1[a]=0;}
for(m=0;m<M;m++){n1[x[m*N+i]]+=W[m];}
for(a=0;a<q;a++){fprintf(fpw,"%lf ",(double)n1[a]/Meff);}
fprintf(fpw,"\n");
}


FILE *fpw2;
fpw2=fopen("stat_align_2p.txt","w");
for(i=0;i<N;i++){
for(j=i+1;j<N;j++){
fprintf(fpw2,"%d %d ",i,j);
for(a=0;a<q*q;a++){n2[a]=0;}
for(m=0;m<M;m++){n2[(x[m*N+i])*q+x[m*N+j]]+=W[m];}
for(a=0;a<q*q;a++){fprintf(fpw2,"%lf ",(double)n2[a]/Meff);}
fprintf(fpw2,"\n");
}}

/*
FILE *fpw3;
fpw3=fopen("stat_align_3p.txt","w");
for(i=0;i<N;i++){
for(j=i+1;j<N;j++){
for(k=j+1;k<N;k++){
fprintf(fpw3,"%d %d %d ",i,j,k);
for(a=0;a<q*q*q;a++){n3[a]=0;}
for(m=0;m<M;m++){n3[(x[m*N+i])*q*q+x[m*N+j]*q+x[m*N+k]]+=W[m];}
for(a=0;a<q*q*q;a++){fprintf(fpw3,"%lf ",(double)n3[a]/Meff);}
fprintf(fpw3,"\n");
}}}
*/

free(x);
free(W);
free(n1);
free(n2);
free(n3);

return 0;
}

