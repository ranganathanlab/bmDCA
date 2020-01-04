#include "mcmc.hpp"

#include <string>

#include "graph.hpp"

void
MCMC::load(potts_model model)
{
  graph.load(model);
}

MCMC::MCMC(size_t N, size_t Q)
  : graph(N, Q)
{
  n = N;
  q = Q;
  srand48(seed);
};

MCMC::MCMC(potts_model params, size_t N, size_t Q)
  : graph(N, Q)
{
  n = N;
  q = Q;
  srand48(seed);
  graph.load(params);
};

void
MCMC::sample(arma::field<arma::Mat<int>>* ptr,
             int reps,
             int M,
             int N,
             int t_wait,
             int delta_t,
             int seed){
#pragma omp parallel
  {
#pragma omp for
    for (int rep = 0; rep < reps; rep++){
      graph.sample_mcmc(&((*ptr).at(rep)), M, t_wait, delta_t, seed + rep);
}
}
}
;

void
MCMC::sample_init(arma::field<arma::Mat<int>>* ptr,
                  int reps,
                  int M,
                  int N,
                  int t_wait,
                  int delta_t,
                  arma::Col<int>* init_ptr){
#pragma omp parallel
  {
#pragma omp for
    for (int rep = 0; rep < reps; rep++){
      graph.sample_mcmc_init(&((*ptr).at(rep)), M, t_wait, delta_t, init_ptr);
}
}
}
;
