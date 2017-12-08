int main(int argc, char const *argv[])
{
	/* code */
	return 0;
}

void read_alignment(char* file_name, char** alignment)
{
	int i,m,f,M,L,q;
	int temp_int;
	FILE* file_handle;
	file_handle = fopen(file_name, "r");
	if(file_handle == NULL)
		{printf("reweighting.c - read_alignment : alignment file %s could not be opened\n",file_name); }
	
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
}