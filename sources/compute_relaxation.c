#include <stdlib.h>
#include <stdio.h>
#include "math.h"


int main(int argc, char *argv[]){

int i,j,k,a,b,m;
int M,N,q;

int RIP,rip;
RIP=atoi(argv[3]);
int DT;
DT=atoi(argv[4]);

int me;
double E_av,E_av2;

//reading sequences...

//printf("lettura..%s\n",argv[1]);
//fflush(0);
FILE *fp;
fp=fopen(argv[1],"r");
fscanf(fp,"%d %d %d\n",&M,&N,&q);


int *x;
x=(int *)malloc(sizeof(int)*N*M);

for(m=0;m<M;m++){
for(i=0;i<N;i++){
fscanf(fp,"%d",&x[m*N+i]);
}
fscanf(fp,"\n");
}



//read parameters
//printf("Read parameters...\n");
//fflush(0);
double E;
int z1,z2,z3,z4;
double *h;
h=(double *)malloc(sizeof(double)*q*N);
double *J;
J=(double *)malloc(sizeof(double)*q*q*N*N);
FILE *fpw;
fpw=fopen(argv[2],"r");

for(i=0;i<N;i++){
for(j=i+1;j<N;j++){
for(a=0;a<q;a++){
for(b=0;b<q;b++){
fscanf(fpw,"J %d %d %d %d %lf\n",&z1,&z2,&z3,&z4,&J[b+a*q+q*q*j+i*N*q*q]);
}}}}
for(i=0;i<N;i++){
for(a=0;a<q;a++){
fscanf(fpw,"h %d %d %lf\n",&z1,&z2,&h[a+i*q]);
}}
fclose(fpw);


//read parameters
//printf("Compute energies...%d\n",M);
//fflush(0);

FILE *fpe;
fpe=fopen(argv[5],"w");

for(me=0;me<M/RIP;me++){
	E_av=0;
	E_av2=0;
	for(rip=0;rip<RIP;rip++){
	m=rip*M/RIP+me;
	E=0;
		for(i=0;i<N;i++){
			E-=h[x[(m)*N+i]+i*q];
			for(j=i+1;j<N;j++)E-=J[x[(m)*N+j]+x[(m)*N+i]*q+q*q*j+i*N*q*q];
			}
	E_av+=E;
	E_av2+=E*E;
	}

fprintf(fpe,"%d %lf %lf\n",me*DT,E_av/RIP,sqrt(E_av2/RIP-pow(E_av/RIP,2))/sqrt(RIP));
}


free(J);
free(h);
free(x);

return 0;
}




