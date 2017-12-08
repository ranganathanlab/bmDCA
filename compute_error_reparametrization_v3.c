#include <stdlib.h>
#include <stdio.h>
#include "math.h"

#define EPSILON 0.00000001 // regularization parameter to estimate deviations (should be something like 1/Meff+1/M)


int main(int argc, char *argv[]){



int i,j,k,a,b,m;
int M,N,q;
double ERROR_MAX,W1,W2;
int z1,z2,z3,z4;
double delta;
double LAMBDA_h,LAMBDA_J;


double rho,beta,den_beta,num_beta,num_rho,den_stat,den_mc,c_mc_av,c_stat_av,rho_1p,num_rho_1p,den_stat_1p,den_mc_1p=0	;
FILE *fp1mc;
fp1mc=fopen(argv[1],"r");
FILE *fp2mc;
fp2mc=fopen(argv[2],"r");
FILE *fp1msa;
fp1msa=fopen(argv[3],"r");
FILE *fp2msa;
fp2msa=fopen(argv[4],"r");
FILE *fperror;
fperror=fopen("errorlist.txt","w");


N=(int)atoi(argv[5]);
q=(int)atoi(argv[6]);
ERROR_MAX=atof(argv[7]);
LAMBDA_h=atof(argv[8]);
LAMBDA_J=atof(argv[9]);
FILE *fp1mcsigma;
fp1mcsigma=fopen(argv[10],"r");
FILE *fp2mcsigma;
fp2mcsigma=fopen(argv[11],"r");

FILE *fpw;
fpw=fopen(argv[12],"r");

double ERROR_MIN_UPDATE;
ERROR_MIN_UPDATE=atof(argv[13]);

double MEFF;
MEFF=atof(argv[14]);

double error_stat_1p,error_stat_2p,error_stat_tot,delta_stat=0;
int count1=0;
int count2=0;

double *n1;
n1=(double *)malloc(sizeof(double)*q*N);
double *n2;
n2=(double *)malloc(sizeof(double)*q*q*N*N);
double *n1mc;
n1mc=(double *)malloc(sizeof(double)*q*N);
double *n2mc;
n2mc=(double *)malloc(sizeof(double)*q*q*N*N);
double *n1mcsigma;
n1mcsigma=(double *)malloc(sizeof(double)*q*N);
double *n2mcsigma;
n2mcsigma=(double *)malloc(sizeof(double)*q*q*N*N);
double *gradh;
gradh=(double *)malloc(sizeof(double)*q*N);
double *gradJ;
gradJ=(double *)malloc(sizeof(double)*q*q*N*N);
double *h;
h=(double *)malloc(sizeof(double)*q*N);
double *J;
J=(double *)malloc(sizeof(double)*q*q*N*N);
double deltamax_1,deltamax_2=0;



/////////////////////////////////////////////////////////////////////////
//reading statistics 1p and 2p
for(i=0;i<N;i++){
	fscanf(fp1msa,"%d ",&z1);
	fscanf(fp1mc,"%d ",&z2);
	fscanf(fp1mcsigma,"%d ",&z2);
		for(a=0;a<q;a++){
		fscanf(fp1msa,"%lf ",&n1[a+q*i]);
		fscanf(fp1mc,"%lf ",&n1mc[a+q*i]);
		fscanf(fp1mcsigma,"%lf ",&n1mcsigma[a+q*i]);
		}
	fscanf(fp1msa,"\n");
	fscanf(fp1mc,"\n");
	fscanf(fp1mcsigma,"\n");
}


for(i=0;i<N;i++){
	for(j=i+1;j<N;j++){
		fscanf(fp2msa,"%d %d ",&z1,&z2);

		fscanf(fp2mc,"%d %d ",&z1,&z2);

		fscanf(fp2mcsigma,"%d %d ",&z1,&z2);

		for(a=0;a<q*q;a++){
			fscanf(fp2msa,"%lf ",&n2[a+q*q*j+i*N*q*q]);
			fscanf(fp2mc,"%lf ",&n2mc[a+q*q*j+i*N*q*q]);
			fscanf(fp2mcsigma,"%lf ",&n2mcsigma[a+q*q*j+i*N*q*q]);
			}
		fscanf(fp2msa,"\n");
		fscanf(fp2mc,"\n");
		fscanf(fp2mcsigma,"\n");
	}
}

//printf("read\n");
//fflush(0);
///////////////////////////////////////////////////////////////////////////////////
//reading parameters


//read
        for(i=0;i<N;i++){
                for(j=i+1;j<N;j++){
                        for(a=0;a<q;a++){
                                for(b=0;b<q;b++){
                                        fscanf(fpw,"J %d %d %d %d %lf\n",&z1,&z2,&z3,&z4,&J[b+a*q+q*q*j+i*N*q*q]);
                                }
                        }
                }
        }



        for(i=0;i<N;i++){
                for(a=0;a<q;a++){
                        fscanf(fpw,"h %d %d %lf\n",&z1,&z2,&h[a+i*q]);
                }
        }

        fclose(fpw);


////////////////////////////////compute error

double error_1p=0;
double error_2p=0;
double error_tot=0;


for(i=0;i<N;i++){

	for(a=0;a<q;a++){
		delta=(n1mc[a+q*i]-n1[a+q*i]+LAMBDA_h*h[a+q*i]);
		delta_stat=(n1mc[a+q*i]-n1[a+q*i])/(pow(n1[a+q*i]*(1-n1[a+q*i])/MEFF+pow(n1mcsigma[a+q*i],2)+EPSILON,0.5)); 		
		error_1p+=pow(delta,2);
		error_stat_1p+=pow(delta_stat,2);
		if (pow(delta,2)>pow(deltamax_1,2))deltamax_1=sqrt(pow(delta,2));
		if(sqrt(pow(delta_stat,2))>ERROR_MIN_UPDATE){
			gradh[a+q*i]=n1[a+q*i]-n1mc[a+q*i]-LAMBDA_h*h[a+q*i];
			count1++;
		}
		fprintf(fperror,"%d %d %lf %lf %lf %lf %lf\n",i,a,n1[a+q*i],n1mc[a+q*i],gradh[a+i*q],(n1[a+q*i]-n1mc[a+q*i])/(pow(n1[a+q*i]*(1-n1[a+q*i])/MEFF+pow(n1mcsigma[a+q*i],2)+EPSILON,0.5)),error_1p);
	}
}

double error_c=0;
double c_mc,c_stat;
int a1,a2;

for(i=0;i<N;i++){

	for(j=i+1;j<N;j++){
		for(a=0;a<q*q;a++){
	
			a1=(int)(a/(q));
			a2=(int)(a-a1*q);			
			delta=-(n2[a+q*q*j+i*N*q*q]-n2mc[a+q*q*j+i*N*q*q]+(n1mc[a1+q*i]-n1[a1+q*i])*n1[a2+q*j]+(n1mc[a2+q*j]-n1[a2+q*j])*n1[a1+q*i]-LAMBDA_J*J[a+q*q*j+i*N*q*q]);
			delta_stat=(n2mc[a+q*q*j+i*N*q*q]-n2[a+q*q*j+i*N*q*q])/(pow(n2[a+q*q*j+i*N*q*q]*(1.0-n2[a+q*q*j+i*N*q*q])/MEFF+pow(n2mcsigma[a+q*q*j+i*N*q*q],2)+EPSILON,0.5));
			
			c_mc=n2mc[a+q*q*j+i*N*q*q]-n1mc[a1+q*i]*n1mc[a2+q*j];
			c_stat=n2[a+q*q*j+i*N*q*q]-n1[a1+q*i]*n1[a2+q*j];
			c_mc_av+=c_mc;
			c_stat_av+=c_stat;
			error_c+=pow(c_mc-c_stat,2);
			error_2p+=pow(delta,2);
			error_stat_2p+=pow(delta_stat,2);

			if (pow(delta,2)>pow(deltamax_2,2)) {deltamax_2=sqrt(pow(delta,2));}
			if(sqrt(pow(delta_stat,2))>ERROR_MIN_UPDATE) {
				gradJ[a+q*q*j+i*N*q*q]=n2[a+q*q*j+i*N*q*q]-n2mc[a+q*q*j+i*N*q*q]+(n1mc[a1+q*i]-n1[a1+q*i])*n1[a2+q*j]+(n1mc[a2+q*j]-n1[a2+q*j])*n1[a1+q*i]-LAMBDA_J*J[a+q*q*j+i*N*q*q];
				count2++;				
				}
		}		
	}
}



c_stat_av/=((N*(N-1)*q*q)/2);
c_mc_av/=((N*(N-1)*q*q)/2);


FILE *fp_corr;
fp_corr=fopen("my_corr.dat","w");

for(i=0;i<N;i++){
	for(j=i+1;j<N;j++){
		for(a=0;a<q*q;a++){
			a1=(int)(a/(q));
			a2=(int)(a-a1*q);
			c_mc=n2mc[a+q*q*j+i*N*q*q]-n1mc[a1+q*i]*n1mc[a2+q*j];
			c_stat=n2[a+q*q*j+i*N*q*q]-n1[a1+q*i]*n1[a2+q*j];
			num_rho+=(c_mc-c_mc_av)*(c_stat-c_stat_av);
			num_beta+=(c_mc)*(c_stat);
			den_stat+=pow(c_stat-c_stat_av,2);			
			den_mc+=pow(c_mc-c_mc_av,2);
			den_beta+=pow(c_stat,2);

fprintf(fp_corr,"%lf %lf\n",c_stat,c_mc);
		}	
	}
}



for(i=0;i<N;i++){
		for(a=0;a<q;a++){
			num_rho_1p+=(n1mc[a+q*i]-1.0/q)*(n1[a+q*i]-1.0/q);			
			den_stat_1p+=pow(n1[a+q*i]-1.0/q,2);
			den_mc_1p+=pow(n1mc[a+q*i]-1.0/q,2);
		}	
}



beta=num_beta/den_beta;
rho=num_rho/sqrt(den_mc*den_stat);
rho_1p=num_rho_1p/sqrt(den_mc_1p*den_stat_1p);



error_1p=sqrt(error_1p/(N*q));
error_2p=sqrt(error_2p/((N*(N-1)*q*q)/2));

error_stat_1p=sqrt(error_stat_1p/(N*q));
error_stat_2p=sqrt(error_stat_2p/(N*(N-1)*q*q)/2);

error_c=sqrt(error_c/(N*(N-1)*q*q)/2);

error_tot=error_1p+error_2p;
error_stat_tot=error_stat_1p+error_stat_2p;



//printf(" error 1p=%lf\n error 2p=%lf\n error tot=%lf EPSILON_h=%lf EPSILON_J=%lf, Errormax1p=%lf Errormax2p=%lf\n",error_1p,error_2p,error_tot,EPSILON_h,EPSILON_J,deltamax_1,deltamax_2);

FILE *fp_error;
fp_error=fopen("error.txt","a");
fprintf(fp_error,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",error_1p,error_2p,error_tot,deltamax_1,deltamax_2,error_stat_1p,error_stat_2p,error_stat_tot,100.0*count1/(double)(N*q),200.0*count2/(double)(N*(N-1)*q*q),error_c,rho,beta,rho_1p);


//printf("error computed\n");
//fflush(0);
///////////////////////////////////////////////////////////////////////////////////////////
// exit or update gradients

if (error_tot<ERROR_MAX){
	printf("converged\n");
	return 0;}
else{
//	printf("update couplings...\n");


	fclose(fpw);
	fpw=fopen("gradient.txt","w");


//update
	for(i=0;i<N;i++){
		for(j=i+1;j<N;j++){
			for(a=0;a<q;a++){
				for(b=0;b<q;b++){
					fprintf(fpw,"J %d %d %d %d %lf\n",i,j,a,b,gradJ[b+a*q+q*q*j+i*N*q*q]);
				}
			}
		}
	}
	for(i=0;i<N;i++){
		for(a=0;a<q;a++){
			fprintf(fpw,"h %d %d %lf\n",i,a,gradh[a+i*q]);
			//printf("h %d %d %lf %lf %lf %lf\n",i,a,deltah[a+i*q],n1[a+q*i],n1mc[a+q*i],log(n1[a+q*i]+LAMBDA_MSA)-log(n1mc[a+q*i]+LAMBDA_MC));
		}
	}
}


//printf("parameters updated\n");
//fflush(0);





free(n1);
free(n2);
free(n1mc);
free(n2mc);
free(n1mcsigma);
free(n2mcsigma);
free(gradh);
free(gradJ);

free(h);
free(J);

return 0;
}

