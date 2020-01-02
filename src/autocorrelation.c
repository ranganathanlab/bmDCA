#include "math.h"
#include <stdio.h>
#include <stdlib.h>

double
min(double x, double y)
{
  if (x < y)
    return x;
  return y;
}
double
max(double x, double y)
{
  if (x < y)
    return y;
  return x;
}

int
main(int argc, char* argv[])
{

  // printf("Computing autocorrelations...\n");
  // fflush(0);

  int i, j, k, a, b, m, m2, id, rip, RIP, rip2;
  double dinf, dinf2;
  int M, N, q;
  FILE* fp;
  fp = fopen(argv[1], "r");
  fscanf(fp, "%d %d %d\n", &M, &N, &q);

  RIP = atoi(argv[2]);

  int DT;
  DT = atoi(argv[3]);

  FILE* fpw;
  fpw = fopen(argv[4], "w");

  int* x;
  x = (int*)malloc(sizeof(int) * N * M);

  double* d;
  d = (double*)malloc(sizeof(double) * (int)(M / RIP));

  double* d2;
  d2 = (double*)malloc(sizeof(double) * (int)(M / RIP));

  int* counter;
  counter = (int*)malloc(sizeof(int) * (int)(M / RIP)); // number of sequences x apart in time

  // printf("M=%d N=%d q=%d\n",M,N,q);
  // fflush(0);

  // reading...
  for (m = 0; m < M; m++) {
    for (i = 0; i < N; i++) {
      fscanf(fp, "%d", &x[m * N + i]);
    }
    fscanf(fp, "\n");
  }

  // compute distances
  /* d[0] = M; */
  /* d2[0] = M; */
  /* counter[0] = M; */
  for (rip = 0; rip < RIP; rip++) {
    for (m = 0; m < M / RIP; m++) {
      for (m2 = m + 1; m2 < M / RIP; m2++) {
        id = 0;
        for (i = 0; i < N; i++) {
          if (x[(rip * M / RIP + m) * N + i] == x[(rip * M / RIP + m2) * N + i])
            id++;
        }
        d[m2 - m] += (double)id / N;
        d2[m2 - m] += (double)id * id / (N * N);
        counter[m2 - m]++;
      }
    }
  }

  dinf = 0;
  dinf2 = 0;

  for (m = 0; m < M / RIP; m++) {
    for (rip = 0; rip < RIP; rip++) {
      for (rip2 = rip + 1; rip2 < RIP; rip2++) {
        id = 0;
        for (i = 0; i < N; i++) {
          if (x[(rip * M / RIP + m) * N + i] == x[(rip2 * M / RIP + m) * N + i])
            id++;
        }
        dinf += (double)id / N;
        dinf2 += (double)id * id / (N * N);
      }
    }
  }

  // output
  for (i = 1; i < M / RIP - 1; i++) {
  /* for (i = 0; i < M / RIP; i++) { */
    fprintf(fpw,
            "%d %lf %lf\n",
            i * DT,
            (double)d[i] / (double)(counter[i]),
            sqrt(1.0 / counter[i]) * sqrt(d2[i] / (double)(counter[i]) -
                                          pow(d[i] / (double)(counter[i]), 2)));
  }

  FILE* fwinf;
  fwinf = fopen("overlap_inf.txt", "w");
  fprintf(fwinf,
          "0 %lf %lf\n",
          (double)2.0 * dinf / (double)(RIP * (RIP - 1) * (M / RIP)),
          sqrt(2.0 / (RIP * (RIP - 1) * (M / RIP))) *
            sqrt(2.0 * dinf2 / (double)(RIP * (RIP - 1) * (M / RIP)) -
                 pow(2.0 * dinf / (double)(RIP * (RIP - 1) * (M / RIP)), 2)));
  fprintf(fwinf,
          "%d %lf %lf\n",
          DT * (M / RIP),
          (double)2.0 * dinf / (double)(RIP * (RIP - 1) * (M / RIP)),
          sqrt(2.0 / (RIP * (RIP - 1) * (M / RIP))) *
            sqrt(2.0 * dinf2 / (double)(RIP * (RIP - 1) * (M / RIP)) -
                 pow(2.0 * dinf / (double)(RIP * (RIP - 1) * (M / RIP)), 2)));

  double sigma_cross, sigma_auto, sigma_check;
  double overlap_auto, overlap_cross, overlap_check;
  double err_cross_auto, err_cross_check, err_check_auto;
  int i_auto, i_check;
  i_check = max((M / RIP) / 10, 1);
  i_auto = 1;

  overlap_cross = (double)2.0 * dinf / (double)(RIP * (RIP - 1) * (M / RIP));
  sigma_cross =
    sqrt(2.0 * dinf2 / (double)(RIP * (RIP - 1) * (M / RIP)) -
         pow(2.0 * dinf / (double)(RIP * (RIP - 1) * (M / RIP)), 2));

  overlap_auto = (double)d[i_auto] / (double)(counter[i_auto]);
  sigma_auto = sqrt(d2[i_auto] / (double)(counter[i_auto]) -
                    pow(d[i_auto] / (double)(counter[i_auto]), 2));

  overlap_check = (double)d[i_check] / (double)(counter[i_check]);
  sigma_check = sqrt(d2[i_check] / (double)(counter[i_check]) -
                     pow(d[i_check] / (double)(counter[i_check]), 2));

  err_cross_auto = sqrt(pow(sigma_cross, 2) + pow(sigma_auto, 2)) / sqrt(RIP);
  err_cross_check = sqrt(pow(sigma_cross, 2) + pow(sigma_check, 2)) / sqrt(RIP);
  err_check_auto = sqrt(pow(sigma_check, 2) + pow(sigma_auto, 2)) / sqrt(RIP);

  FILE* fp_ergo;
  fp_ergo = fopen("ergo.txt", "a");

  fprintf(fp_ergo,
          "%lf %lf %lf %lf %lf %lf %lf %lf %f\n",
          overlap_auto,
          overlap_check,
          overlap_cross,
          sigma_auto,
          sigma_check,
          sigma_cross,
          err_cross_auto,
          err_cross_check,
          err_check_auto);

  free(x);
  free(counter);
  free(d2);
  free(d);

  return 0;
}
