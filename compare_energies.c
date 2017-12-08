#include <stdlib.h>
#include <stdio.h>
#include "math.h"


int main(int argc, char *argv[]){

int i,m;
int M,N,q;


double delta2,rho,rho_num,den1,den2,e1_av,e2_av=0;


//reading sequences...
FILE *fp;
fp=fopen(argv[1],"r");
FILE *fp2;
fp2=fopen(argv[2],"r");
M=atoi(argv[3]);

FILE *fp_out;
fp_out=fopen(argv[4],"a");

double *e1;
e1=(double *)malloc(sizeof(double)*M);
double *e2;
e2=(double *)malloc(sizeof(double)*M);

for (i=0;i<M;i++){
fscanf(fp,"%lf ",&e1[i]);

fscanf(fp2,"%lf ",&e2[i]);
e1_av+=e1[i];
e2_av+=e2[i];
delta2+=pow(e1[i]-e2[i],2);
}

delta2=sqrt(delta2/M);
e1_av/=M;
e2_av/=M;



for (i=0;i<M;i++){
rho_num+=(e1[i]-e1_av)*(e2[i]-e2_av);
den1+=pow(e1[i]-e1_av,2);
den2+=pow(e2[i]-e2_av,2);
}

rho=rho_num/sqrt(den1*den2);

fprintf(fp_out,"%lf %lf\n",rho,delta2);


free(e1);
free(e2);

return 0;
}
