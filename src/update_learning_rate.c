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

int
Theta(double x)
{
  if (x > 0)
    return 1;
  return 0;
}

int
Delta(double x)
{
  if (x == 0)
    return 1;
  return 0;
}

double
Max(double a, double b)
{
  if (a > b)
    return a;
  return b;
}

double
Min(double a, double b)
{
  if (a < b)
    return a;
  return b;
}

int
main(int argc, char* argv[])
{

  int i, j, k, a, b, m;
  int M, N, q;
  double ADAPT_UP, ADAPT_DOWN;
  int z1, z2, z3, z4;

  N = (int)atoi(argv[1]); // n samples
  q = (int)atoi(argv[2]); // alphabet size
  char* learning_file;
  learning_file = argv[3]; // learning_rate.txt
  char* gradient_file;
  gradient_file = argv[4]; // gradient.txt
  char* gradient_file_old;
  gradient_file_old = argv[5]; // gradient_old.txt

  ADAPT_UP = atof(argv[6]);
  ADAPT_DOWN = atof(argv[7]);

  double MIN_STEPh;
  MIN_STEPh = atof(argv[8]);

  double MAX_STEPh;
  MAX_STEPh = atof(argv[9]);

  double MIN_STEPj;
  MIN_STEPj = atof(argv[10]);

  double MAX_STEPj;
  MAX_STEPj = atof(argv[11]);
  // TOLERANCE=atof(argv[8]);

  ///////////read learning rates and gradients

  double* eps_h;
  eps_h = (double*)malloc(sizeof(double) * q * N);
  double* eps_J;
  eps_J = (double*)malloc(sizeof(double) * q * q * N * N);
  FILE* fpl;
  fpl = fopen(learning_file, "r");

  // FILE *fpl2;
  // fpl2=fopen("history_sign.txt","a");

  double* gradh;
  gradh = (double*)malloc(sizeof(double) * q * N);
  double* gradJ;
  gradJ = (double*)malloc(sizeof(double) * q * q * N * N);
  FILE* fpg;
  fpg = fopen(gradient_file, "r");

  double* gradh_old;
  gradh_old = (double*)malloc(sizeof(double) * q * N);
  double* gradJ_old;
  gradJ_old = (double*)malloc(sizeof(double) * q * q * N * N);
  FILE* fpg_old;
  fpg_old = fopen(gradient_file_old, "r");

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
          fscanf(fpg,
                 "J %d %d %d %d %lf\n",
                 &z1,
                 &z2,
                 &z3,
                 &z4,
                 &gradJ[b + a * q + q * q * j + i * N * q * q]);
          fscanf(fpg_old,
                 "J %d %d %d %d %lf\n",
                 &z1,
                 &z2,
                 &z3,
                 &z4,
                 &gradJ_old[b + a * q + q * q * j + i * N * q * q]);
        }
      }
    }
  }

  for (i = 0; i < N; i++) {
    for (a = 0; a < q; a++) {
      fscanf(fpl, "h %d %d %lf\n", &z1, &z2, &eps_h[a + i * q]);
      fscanf(fpg, "h %d %d %lf\n", &z1, &z2, &gradh[a + i * q]);
      fscanf(fpg_old, "h %d %d %lf\n", &z1, &z2, &gradh_old[a + i * q]);
    }
  }
  // printf("%lf %lf\n",gradJ[q*q],gradJ_old[q*q]);

  fclose(fpl);
  fclose(fpg);
  fclose(fpg_old);

  // write
  fpl = fopen(learning_file, "w");
  double alfa;

  for (i = 0; i < N; i++) {
    for (j = i + 1; j < N; j++) {
      for (a = 0; a < q; a++) {
        for (b = 0; b < q; b++) {
          alfa = Theta(gradJ[b + a * q + q * q * j + i * N * q * q] *
                       gradJ_old[b + a * q + q * q * j + i * N * q * q]) *
                   ADAPT_UP +
                 Theta(-gradJ[b + a * q + q * q * j + i * N * q * q] *
                       gradJ_old[b + a * q + q * q * j + i * N * q * q]) *
                   ADAPT_DOWN +
                 Delta(gradJ[b + a * q + q * q * j + i * N * q * q] *
                       gradJ_old[b + a * q + q * q * j + i * N * q * q]);
          // fprintf(fpl2,"%d
          // ",Theta(gradJ[b+a*q+q*q*j+i*N*q*q]*gradJ_old[b+a*q+q*q*j+i*N*q*q])-Theta(-gradJ[b+a*q+q*q*j+i*N*q*q]*gradJ_old[b+a*q+q*q*j+i*N*q*q]));

          fprintf(
            fpl,
            "J %d %d %d %d %lf\n",
            i,
            j,
            a,
            b,
            Min(MAX_STEPj,
                Max(MIN_STEPj,
                    alfa * eps_J[b + a * q + q * q * j + i * N * q * q])));
        }
      }
    }
  }

  for (i = 0; i < N; i++) {
    for (a = 0; a < q; a++) {
      alfa = Theta(gradh[a + i * q] * gradh_old[a + i * q]) * ADAPT_UP +
             Theta(-gradh[a + i * q] * gradh_old[a + i * q]) * ADAPT_DOWN +
             Delta(gradh[a + i * q] * gradh_old[a + i * q]);
      // fprintf(fpl2,"%d
      // ",Theta(gradh[a+i*q]*gradh_old[a+i*q])-Theta(-gradh[a+i*q]*gradh_old[a+i*q]));
      fprintf(fpl,
              "h %d %d %lf\n",
              i,
              a,
              Min(MAX_STEPh, Max(MIN_STEPh, alfa * eps_h[a + i * q])));
    }
  }

  // fprintf(fpl2,"\n");
  // fclose(fpl2);
  fclose(fpl);

  ///////////////////////

  free(gradh_old);
  free(gradJ_old);
  free(gradh);
  free(gradJ);
  free(eps_h);
  free(eps_J);

  return 0;
}
