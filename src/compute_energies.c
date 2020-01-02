#include "math.h"
#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char* argv[])
{

  int i, j, k, a, b, m;
  int M, N, q;

  // reading sequences...

  FILE* fp;
  fp = fopen(argv[1], "r");
  fscanf(fp, "%d %d %d\n", &M, &N, &q);

  int* x;
  x = (int*)malloc(sizeof(int) * N * M);

  for (m = 0; m < M; m++) {
    for (i = 0; i < N; i++) {
      fscanf(fp, "%d", &x[m * N + i]);
    }
    fscanf(fp, "\n");
  }

  // read parameters
  double E;
  int z1, z2, z3, z4;
  double* h;
  h = (double*)malloc(sizeof(double) * q * N);
  double* J;
  J = (double*)malloc(sizeof(double) * q * q * N * N);
  FILE* fpw;
  fpw = fopen(argv[2], "r");

  for (i = 0; i < N; i++) {
    for (j = i + 1; j < N; j++) {
      for (a = 0; a < q; a++) {
        for (b = 0; b < q; b++) {
          fscanf(fpw,
                 "J %d %d %d %d %lf\n",
                 &z1,
                 &z2,
                 &z3,
                 &z4,
                 &J[b + a * q + q * q * j + i * N * q * q]);
        }
      }
    }
  }
  for (i = 0; i < N; i++) {
    for (a = 0; a < q; a++) {
      fscanf(fpw, "h %d %d %lf\n", &z1, &z2, &h[a + i * q]);
    }
  }
  fclose(fpw);

  FILE* fpe;
  fpe = fopen(argv[3], "w");
  for (m = 0; m < M; m++) {
    E = 0;
    for (i = 0; i < N; i++) {
      E -= h[x[(m)*N + i] + i * q];
      for (j = i + 1; j < N; j++)
        E -= J[x[(m)*N + j] + x[(m)*N + i] * q + q * q * j + i * N * q * q];
    }
    fprintf(fpe, "%lf\n", E);
  }

  free(J);
  free(h);
  free(x);

  return 0;
}
