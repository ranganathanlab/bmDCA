//v3= ho inserito il check sui rapporti tra le z
//v4 ho inserito i pesi con l'ipr

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


FILE *fp_importance;
fp_importance=fopen("coherence_importance.txt","a");

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

double *n1;
n1=(double *)malloc(sizeof(double)*q*RIP);

double *n2;
n2=(double *)malloc(sizeof(double)*q*q*RIP);


double *n1av;
n1av=(double *)malloc(sizeof(double)*q);

double *n2av;
n2av=(double *)malloc(sizeof(double)*q*q);

double *n1squared;
n1squared=(double *)malloc(sizeof(double)*q);

double *n2squared;
n2squared=(double *)malloc(sizeof(double)*q*q);



//printf("M=%d N=%d q=%d MTOT=%d\n",M,N,q,MTOT);
//fflush(0);


//reading...
for(m=0;m<MTOT;m++){
for(i=0;i<N;i++){
fscanf(fp,"%d",&x[m*N+i]);
}
fscanf(fp,"\n");
}

//////////////read parameters
  double *h;
        h=(double *)malloc(sizeof(double)*q*N);
        double *J;
        J=(double *)malloc(sizeof(double)*q*q*N*N);
        FILE *fp_par;
        fp_par=fopen(argv[3],"r");
 
 double *h_ref;
        h_ref=(double *)malloc(sizeof(double)*q*N);
        double *J_ref;
        J_ref=(double *)malloc(sizeof(double)*q*q*N*N);
        FILE *fp_par_ref;
        fp_par_ref=fopen(argv[4],"r");

int z1,z2,z3,z4;
//read
        for(i=0;i<N;i++){
                for(j=i+1;j<N;j++){
                        for(a=0;a<q;a++){
                                for(b=0;b<q;b++){
                                        fscanf(fp_par,"J %d %d %d %d %lf\n",&z1,&z2,&z3,&z4,&J[b+a*q+q*q*j+i*N*q*q]);
					fscanf(fp_par_ref,"J %d %d %d %d %lf\n",&z1,&z2,&z3,&z4,&J_ref[b+a*q+q*q*j+i*N*q*q]);
                                }
                        }
                }
        }



        for(i=0;i<N;i++){
                for(a=0;a<q;a++){
                        fscanf(fp_par,"h %d %d %lf\n",&z1,&z2,&h[a+i*q]);
			fscanf(fp_par_ref,"h %d %d %lf\n",&z1,&z2,&h_ref[a+i*q]);
                }
        }

        fclose(fp_par);

//////compute energy
double *dE;
double *p;
double *dE_av;
double *Z,*Z_inv,Z_tot,Z_inv_tot;
dE=(double *)malloc(sizeof(double)*MTOT);
dE_av=(double *)malloc(sizeof(double)*RIP);
Z=(double *)malloc(sizeof(double)*RIP);
Z_inv=(double *)malloc(sizeof(double)*RIP);
p=(double *)malloc(sizeof(double)*MTOT);
double *sum;
sum=(double *)malloc(sizeof(double)*RIP);
double *w;
w=(double *)malloc(sizeof(double)*RIP);
Z_tot=0;
Z_inv_tot=0;
double W=0;
double sumw=0;
double dE_av_tot=0;

for(rip=0;rip<RIP;rip++){	
	dE_av[rip]=0;
	for(m=0;m<M;m++){
		dE[m+M*rip]=0;
		for(i=0;i<N;i++){
			dE[m+M*rip]+=h[x[(rip*M+m)*N+i]+i*q]-h_ref[x[(rip*M+m)*N+i]+i*q];
			for(j=i+1;j<N;j++)dE[m+M*rip]+=J[x[(rip*M+m)*N+j]+x[(rip*M+m)*N+i]*q+q*q*j+i*N*q*q]-J_ref[x[(rip*M+m)*N+j]+x[(rip*M+m)*N+i]*q+q*q*j+i*N*q*q];
			}
		dE_av[rip]+=dE[m+M*rip];
		}
	dE_av[rip]/=M;
	dE_av_tot+=dE_av[rip];
	}

for(rip=0;rip<RIP;rip++){
	Z[rip]=0;
	Z_inv[rip]=0;
	sum[rip]=0;
	for(m=0;m<M;m++){
		p[rip*M+m]=exp((dE[rip*M+m]-dE_av[rip]));
		Z[rip]+=p[rip*M+m];
		Z_inv[rip]+=1/p[rip*M+m];
	}
	//printf("%d %lf %lf **",rip,Z[rip],dE_av[rip]);
	for(m=0;m<M;m++){
		p[rip*M+m]/=(Z[rip]);
		sum[rip]+=pow(p[rip*M+m],2);
		//printf("%d %d %lf %lf %lf\n",rip,m,p[rip*M+m],(dE[rip*M+m]-dE_av[rip]),dE[rip*M+m]);
	}
	Z_tot+=Z[rip];
	Z_inv_tot+=Z_inv[rip];
	w[rip]=1/sum[rip];
	W+=w[rip];
}


for(rip=0;rip<RIP;rip++){
w[rip]=w[rip]/W;
sumw+=pow(w[rip],2);
printf("%d %lf %lf %lf\n",rip,w[rip],Z[rip]/Z_inv[rip],dE_av[rip]);
}

printf("ZTOT ratio %lf, sumw %lf,dE_av_tot=%lf\n",Z_tot/Z_inv_tot,1/sumw,dE_av_tot);

fprintf(fp_importance,"%lf %lf %lf\n",Z_tot/Z_inv_tot,1/sumw,dE_av_tot);
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
n1[rip*q+x[(rip*M+m)*N+i]]+=p[m+M*rip];
}
for(a=0;a<q;a++){
for(rip=0;rip<RIP;rip++){
n1av[a]+=w[rip]*n1[rip*q+a];
n1squared[a]+=w[rip]*pow(n1[rip*q+a],2);
}
//printf("%d %d %lf %lf %lf %lf %d %d %d\n",i,a,(double)n1squared[a]/(M*M*RIP),pow((double)n1av[a]/(M*RIP),2),(double)n1av[a]/(M*RIP),sqrt(((double)n1squared[a]/(M*M*RIP)-pow((double)n1av[a]/(M*RIP),2))/sqrt(RIP)),n1[0*q+a],n1[1*q+a],n1[2*q+a]);
fprintf(fpw,"%lf ",(double)n1av[a]);
fprintf(fpwsigma,"%lf ",max(sqrt(((double)n1squared[a]-pow((double)n1av[a],2))*sqrt(sumw)),0));}
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
for(rip=0;rip<RIP;rip++)for(m=0;m<M;m++){n2[rip*q*q+(x[(rip*M+m)*N+i])*q+x[(rip*M+m)*N+j]]+=p[m+rip*M];}

for(a=0;a<q*q;a++){
for(rip=0;rip<RIP;rip++){
n2av[a]+=w[rip]*n2[rip*q*q+a];
n2squared[a]+=w[rip]*pow(n2[rip*q*q+a],2);
}

//if(j==9)printf("%d %d %d %lf %lf %lf %lf %d %d %d\n",j,i,a,(double)n2squared[a]/(M*M*RIP),pow((double)n2av[a]/(M*RIP),2),(double)n2av[a]/(M*RIP),sqrt(((double)n2squared[a]/(M*M*RIP)-pow((double)n2av[a]/(M*RIP),2))/sqrt(RIP)),n2[0*q*q+a],n2[1*q*q+a],n2[2*q*q+a]);
fprintf(fpw2,"%lf ",(double)n2av[a]);
fprintf(fpw2sigma,"%lf ",max(sqrt(((double)n2squared[a]-pow((double)n2av[a],2))*sqrt(sumw)),0));}




fprintf(fpw2,"\n");
fprintf(fpw2sigma,"\n");
}}

free(w);
free(Z_inv);
free(Z);
free(sum);
free(J);
free(J_ref);
free(h);
free(h_ref);
free(p);
free(dE);
free(dE_av);
free(x);
free(n1);
free(n2);
free(n1av);
free(n2av);
free(n1squared);
free(n2squared);


return 0;
}

