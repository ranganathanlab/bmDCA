#include <stdlib.h>
#include <stdio.h>
#include "math.h"



int main(int argc, char *argv[]){

int i,j,a,b,m;
int N,q;
double xJ,xh;

N=(int)atoi(argv[1]);
q=(int)atoi(argv[2]);
char *learning_file;
learning_file=argv[3];
xJ=atof(argv[4]);
xh=atof(argv[5]);


///////////read learning rates and gradie


// write
	FILE *fpl;
	fpl=fopen(learning_file,"w");
	

     for(i=0;i<N;i++){
                for(j=i+1;j<N;j++){
                        for(a=0;a<q;a++){
                                for(b=0;b<q;b++){
                                        fprintf(fpl,"J %d %d %d %d %lf\n",i,j,a,b,xJ);
                                }
                        }
                }
        }



        for(i=0;i<N;i++){
                for(a=0;a<q;a++){
                        fprintf(fpl,"h %d %d %lf\n",i,a,xh);
                }
        }

        fclose(fpl);

///////////////////////

return 0;
}


