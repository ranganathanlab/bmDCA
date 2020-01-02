#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char* argv[])
{

  int i, j, a, b, m, z1;
  int N, q;
  double pc;
  char* output_file;
  char* stat_file;
  double* n1;

  N = (int)atoi(argv[1]); // number of positions
  q = (int)atoi(argv[2]); // number of amino acids
  output_file = argv[3]; // parameters_temp.txt
  stat_file = argv[4]; // read in 1p statistics
  pc = (double)atof(argv[5]); // 1 / M_eff (i.e. 1 / effective # of sequences)

  n1 = (double*)malloc(sizeof(double) * N * q);

  // read
  FILE* fp1msa;
  fp1msa = fopen(stat_file, "r");

  for (i = 0; i < N; i++) {
    fscanf(fp1msa, "%d ", &z1);
    for (a = 0; a < q; a++) {

      fscanf(fp1msa, "%lf", &n1[a + q * i]);
    }
    fscanf(fp1msa, "\n");
  }

  // write
  FILE* fpl;
  fpl = fopen(output_file, "w");

  for (i = 0; i < N; i++) {
    for (j = i + 1; j < N; j++) {
      for (a = 0; a < q; a++) {
        for (b = 0; b < q; b++) {
          fprintf(fpl, "J %d %d %d %d 0\n", i, j, a, b);
        }
      }
    }
  }

  double avg;
  for (i = 0; i < N; i++) {
    avg = 0;
    for (a = 0; a < q; a++) {
      avg += log((1 - pc) * n1[a + q * i] + pc / q);
    }
    for (a = 0; a < q; a++) {
      fprintf(fpl,
              "h %d %d %lf\n",
              i,
              a,
              log((1 - pc) * n1[a + q * i] + pc / q) - avg / q);
    }
  }

  fclose(fpl);
  fclose(fp1msa);
  free(n1);

  printf("pseudocount: %lf\n", pc);
  fflush(0);

  return 0;
}
