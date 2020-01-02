#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MEFF                                                                   \
  1000.0 // regularization parameter to estimate deviations (should be something
         // like 1/Meff+1/M)
#define MMC                                                                    \
  10000.0 // regularization parameter to estimate deviations (should be
          // something like 1/Meff+1/M)
#define EPSILON                                                                \
  0.0001 // regularization parameter to estimate deviations (should be something
         // like 1/Meff+1/M)
#define ERROR_MIN_UPDATE 0

int
main(int argc, char* argv[])
{

  int i, j, k, a, b, m;
  int M, N, q;
  double ERROR_MAX, EPSILON_J, EPSILON_h, W1, W2;
  int z1, z2, z3, z4;
  double delta;
  double LAMBDA_h, LAMBDA_J;

  N = (int)atoi(argv[1]);
  q = (int)atoi(argv[2]);
  char* learning_file = argv[3];
  char* parameter_file = argv[4];
  char* gradient_file = argv[5];
  char* output_file = argv[6];
  char* stat1p_file = argv[7];

  double* gradh;
  gradh = (double*)malloc(sizeof(double) * q * N);
  double* gradJ;
  gradJ = (double*)malloc(sizeof(double) * q * q * N * N);

  ///////////read learning rates

  double* n1;
  n1 = (double*)malloc(sizeof(double) * q * N);
  double* Dh;
  Dh = (double*)malloc(sizeof(double) * q * N);
  double* eps_h;
  eps_h = (double*)malloc(sizeof(double) * q * N);
  double* eps_J;
  eps_J = (double*)malloc(sizeof(double) * q * q * N * N);
  FILE* fpl;
  fpl = fopen(learning_file, "r");

  // read
  for (i = 0; i < N; i++) {
    for (j = i + 1; j < N; j++) {
      for (a = 0; a < q; a++) {
        for (b = 0; b < q; b++) {
          fscanf(fpl,
                 "J %d %d %d %d %lf\n",
                 &z1,
                 &z2,
                 &z3,
                 &z4,
                 &eps_J[b + a * q + q * q * j + i * N * q * q]);
        }
      }
    }
  }

  for (i = 0; i < N; i++) {
    for (a = 0; a < q; a++) {
      fscanf(fpl, "h %d %d %lf\n", &z1, &z2, &eps_h[a + i * q]);
    }
  }

  fclose(fpl);

  /////////////////read parameters

  double* h;
  h = (double*)malloc(sizeof(double) * q * N);
  double* J;
  J = (double*)malloc(sizeof(double) * q * q * N * N);
  FILE* fpw;
  fpw = fopen(parameter_file, "r");

  // read
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

  /////////////////read gradient

  fpw = fopen(gradient_file, "r");

  // read
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
                 &gradJ[b + a * q + q * q * j + i * N * q * q]);
        }
      }
    }
  }

  for (i = 0; i < N; i++) {
    for (a = 0; a < q; a++) {
      fscanf(fpw, "h %d %d %lf\n", &z1, &z2, &gradh[a + i * q]);
    }
  }

  fclose(fpw);

  ///////////////read marginals

  fpw = fopen(stat1p_file, "r");
  // reading
  for (i = 0; i < N; i++) {
    fscanf(fpw, "%d ", &z1);
    for (a = 0; a < q; a++) {
      fscanf(fpw, "%lf ", &n1[a + q * i]);
    }
    fscanf(fpw, "\n");
  }
  fclose(fpw);
  ////////////////////////////////update

  fpw = fopen(output_file, "w");

  // update
  for (i = 0; i < N; i++) {
    for (j = i + 1; j < N; j++) {
      for (a = 0; a < q; a++) {
        for (b = 0; b < q; b++) {
          fprintf(fpw,
                  "J %d %d %d %d %lf\n",
                  i,
                  j,
                  a,
                  b,
                  J[b + a * q + q * q * j + i * N * q * q] +
                    eps_J[b + a * q + q * q * j + i * N * q * q] *
                      gradJ[b + a * q + q * q * j + i * N * q * q]);
        }
      }
    }
  }

  for (i = 0; i < N; i++) {
    for (a = 0; a < q; a++) {
      for (j = 0; j < N; j++) {
        if (i < j) {
          for (b = 0; b < q; b++) {
            Dh[a + i * q] += -n1[b + j * q] *
                             eps_J[b + a * q + q * q * j + i * N * q * q] *
                             gradJ[b + a * q + q * q * j + i * N * q * q];
          }
        }
        if (i > j) {
          for (b = 0; b < q; b++) {
            Dh[a + i * q] += -n1[b + j * q] *
                             eps_J[a + b * q + q * q * i + j * N * q * q] *
                             gradJ[a + b * q + q * q * i + j * N * q * q];
          }
        }
      }
    }
  }

  for (i = 0; i < N; i++) {
    for (a = 0; a < q; a++) {
      fprintf(fpw,
              "h %d %d %lf\n",
              i,
              a,
              h[a + i * q] + eps_h[a + i * q] * gradh[a + i * q] +
                Dh[a + i * q]);
    }
  }

  // printf("parameters updated\n");
  // fflush(0);

  free(n1);
  free(Dh);
  free(gradh);
  free(gradJ);
  free(eps_h);
  free(eps_J);
  free(h);
  free(J);

  return 0;
}
