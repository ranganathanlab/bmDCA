#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void write_weights(char* file_name, double* weights, int M);
char** read_alignment(char* file_name, int* M, int* L);
double* compute_weights(char** alignment, int M, int L, double threshold);

int main(int argc, char const *argv[])
{
	char in_file[500] = "";
	char out_file[500] = "";
	char** alignment;
	double* weights;
	int M,L;
	double theta = 0.8;

	strcpy(in_file,argv[1]);
	strcpy(out_file, argv[2]);
	printf("	Reading input alignment ...");
	alignment = read_alignment(in_file, &M, &L);
	printf("	Done.\n");
	printf("	Computing weights ... ");
	fflush(stdout);
	weights = compute_weights(alignment, M, L, theta);
	printf("	Done\n");
	printf("	Writing output ... ");
	write_weights(out_file, weights,M);
	printf("	Done\n");



	return 0;
}

void write_weights(char* file_name, double* weights, int M)
{
	int m;
	FILE* handle;
	// printf("write_weights - M=%d -- out file = %s\n",M,file_name);
	handle = fopen(file_name,"w");
	for (m = 0; m < M; ++m)
	{	
		// printf("m = %d/%d\n",m,M);
		fprintf(handle,"%f\n",weights[m]);
	}
	fclose(handle);
}	

char** read_alignment(char* file_name, int* pM, int* pL)
{
	int i,m,f,q,temp_int,M,L;
	char** alignment = NULL;
	FILE* file_handle;

	// Opening file
	file_handle = fopen(file_name, "r");
	if(file_handle == NULL)
		{printf("reweighting.c - read_alignment : alignment file %s could not be opened\n",file_name); }

	// Reading header
	f = fscanf(file_handle,"%d %d %d\n", pM, pL,&q);
	if(f!=3)
		{printf("reweighting.c - read_alignment : There may be a problem with header of alignment file %s.\n",file_name);}
	M = *pM;
	L = *pL;
	// Allocating alignment
	alignment = malloc(M*sizeof(*alignment));
	for (m = 0; m < M; ++m)
	{
		alignment[m] = malloc(L*sizeof(**alignment));	
	}	
	// Reading alignment
	for (m = 0; m < M; ++m)
	{
		for (i = 0; i < L; ++i)
		{	
			if(feof(file_handle))
			{
				break;
			}

			f = fscanf(file_handle,"%d ",&temp_int);
			if(f!=1)
				{printf("reweighting.c - read_alignment : WARNING fscanf could not read one integer here (m=%d, i=%d)\n",m,i);}
			alignment[m][i] = (char) temp_int-1;
		}
	}
	if( (!feof(file_handle)) & (fgetc(file_handle)!=EOF) & (!feof(file_handle)))
	{	
		printf("%c %c %c\n",fgetc(file_handle),fgetc(file_handle),fgetc(file_handle));
		printf("reweighting.c - read_alignment : WARNING alignment file %s may not have the correct number of elements (EOF not reached)\n", file_name); 
		
	}
	
	fclose(file_handle);
	return alignment;
}




double* compute_weights(char** alignment, int M, int L, double threshold)
{
	int m1,m2,i,comp;
	double id;
	double* weights;
	
	weights = malloc(M*sizeof(*weights));
	for (m1 = 0; m1 < M; ++m1)
	{
		weights[m1] = 0;
	}
	comp = 0;
	printf("\n");
	for (m1 = 0; m1 < M; ++m1)
	{
		if((100*m1)/M > comp)
		{	
			fflush(stdout);
			printf("	Completion : %d%%												\r",comp);
			fflush(stdout);
			comp+=5;
		}
		weights[m1] ++;
		for (m2 = m1+1; m2 < M; ++m2)
		{	
			id=0;
			for (i = 0; i < L; ++i)
			{
				if(alignment[m1][i]==alignment[m2][i])
				{
					id++;
				}
			}
			if(id>threshold*L)
			{

				weights[m1] ++;
				weights[m2] ++;
			}
		}
	}
	printf("\n");
	for (m1 = 0; m1 < M; ++m1)
	{
		weights[m1] = 1./weights[m1];
	}
	return weights;
}
