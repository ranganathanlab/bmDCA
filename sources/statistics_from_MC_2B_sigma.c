#include <stdlib.h>
#include <stdio.h>
#include "math.h"



double max(double a,double b){
if (a>b) return a;
return b;
}

int main(int argc, char *argv[]){

//printf("computing statistics MC...\n");
//fflush(0);


int i,j,k,a,b,m,rip;
int M,N,q,MTOT;
FILE *fp;
fp=fopen(argv[1],"r");
fscanf(fp,"%d %d %d\n",&MTOT,&N,&q);

int RIP;
RIP=atoi(argv[2]);

M=MTOT/RIP;

int *x;
x=(int *)malloc(sizeof(int)*N*MTOT);

int *n1;
n1=(int *)malloc(sizeof(int)*q*RIP);

int *n2;
n2=(int *)malloc(sizeof(int)*q*q*RIP);


int *n1av;
n1av=(int *)malloc(sizeof(int)*q);

int *n2av;
n2av=(int *)malloc(sizeof(int)*q*q);

int *n1squared;
n1squared=(int *)malloc(sizeof(int)*q);

int *n2squared;
n2squared=(int *)malloc(sizeof(int)*q*q);



//printf("M=%d N=%d q=%d MTOT=%d\n",M,N,q,MTOT);
//fflush(0);


//reading...
for(m=0;m<MTOT;m++){
for(i=0;i<N;i++){
fscanf(fp,"%d",&x[m*N+i]);
}
fscanf(fp,"\n");
}


//output
FILE *fpw;
fpw=fopen("stat_MC_1p.txt","w");
FILE *fpwsigma;
fpwsigma=fopen("stat_MC_1p_sigma.txt","w");
for(i=0;i<N;i++){
fprintf(fpw,"%d ",i);
fprintf(fpwsigma,"%d ",i);
for(a=0;a<q*RIP;a++){n1[a]=0;}
for(a=0;a<q;a++){n1av[a]=0;
n1squared[a]=0;}
for(rip=0;rip<RIP;rip++)for(m=0;m<M;m++){
n1[rip*q+x[(rip*M+m)*N+i]]++;
}
for(a=0;a<q;a++){
for(rip=0;rip<RIP;rip++){
n1av[a]+=n1[rip*q+a];
n1squared[a]+=pow(n1[rip*q+a],2);
}
//printf("%d %d %lf %lf %lf %lf %d %d %d\n",i,a,(double)n1squared[a]/(M*M*RIP),pow((double)n1av[a]/(M*RIP),2),(double)n1av[a]/(M*RIP),sqrt(((double)n1squared[a]/(M*M*RIP)-pow((double)n1av[a]/(M*RIP),2))/sqrt(RIP)),n1[0*q+a],n1[1*q+a],n1[2*q+a]);
fprintf(fpw,"%lf ",(double)n1av[a]/(M*RIP));
fprintf(fpwsigma,"%lf ",max(sqrt(((double)n1squared[a]/(M*M*RIP)-pow((double)n1av[a]/(M*RIP),2))/sqrt(RIP)),0));}
fprintf(fpw,"\n");
fprintf(fpwsigma,"\n");
}


FILE *fpw2;
fpw2=fopen("stat_MC_2p.txt","w");
FILE *fpw2sigma;
fpw2sigma=fopen("stat_MC_2p_sigma.txt","w");

for(i=0;i<N;i++){
//printf("%d ",i);
for(j=i+1;j<N;j++){
fprintf(fpw2,"%d %d ",i,j);
fprintf(fpw2sigma,"%d %d ",i,j);
for(a=0;a<q*q*RIP;a++){n2[a]=0;}
for(a=0;a<q*q;a++){n2av[a]=0;
n2squared[a]=0;}
for(rip=0;rip<RIP;rip++)for(m=0;m<M;m++){n2[rip*q*q+(x[(rip*M+m)*N+i])*q+x[(rip*M+m)*N+j]]++;}

for(a=0;a<q*q;a++){
for(rip=0;rip<RIP;rip++){
n2av[a]+=n2[rip*q*q+a];
n2squared[a]+=pow(n2[rip*q*q+a],2);
}
//if(j==9)printf("%d %d %d %lf %lf %lf %lf %d %d %d\n",j,i,a,(double)n2squared[a]/(M*M*RIP),pow((double)n2av[a]/(M*RIP),2),(double)n2av[a]/(M*RIP),sqrt(((double)n2squared[a]/(M*M*RIP)-pow((double)n2av[a]/(M*RIP),2))/sqrt(RIP)),n2[0*q*q+a],n2[1*q*q+a],n2[2*q*q+a]);
fprintf(fpw2,"%lf ",(double)n2av[a]/(M*RIP));
fprintf(fpw2sigma,"%lf ",max(sqrt(((double)n2squared[a]/(M*M*RIP)-pow((double)n2av[a]/(M*RIP),2))/sqrt(RIP)),0));}
fprintf(fpw2,"\n");
fprintf(fpw2sigma,"\n");
}}



free(x);
free(n1);
free(n2);
free(n1av);
free(n2av);
free(n1squared);
free(n2squared);


return 0;
}

