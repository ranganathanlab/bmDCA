#include "mcmc.h"

#include "graph.hpp"

// void MCMC::initializeGraph(potts_model model) {
// };

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
  rip = 10;
  m = 10;
  mtot = 10;
  seed = 1;
  mc_iters0 = 1000000;
  mc_iters = 100000;
  mc_rip = 1;
  gauge0 = false;
  // init = false;
  suffix = "";

  srand48(seed);
};

MCMC::MCMC(potts_model params, size_t N, size_t Q)
  : graph(N, Q)
{
  rip = 10;
  n = N;
  q = Q;
  m = 10;
  mtot = 10;
  seed = 1;
  mc_iters0 = 1000000;
  mc_iters = 100000;
  mc_rip = 1;
  gauge0 = false;
  // init = false;
  suffix = "";

  srand48(seed);
  mtot = m * mc_rip;

  graph.load(params);
};

void
MCMC::sample(arma::field<arma::Mat<int>>* ptr,
             int reps,
             int M,
             int N,
             int mc_iters0,
             int mc_iters,
             int seed)
{
#pragma omp parallel
{
  // std::cout << omp_get_num_threads() << std::endl;
#pragma omp for
  for (int rep = 0; rep < reps; rep++) {
    graph.sample_mcmc(&((*ptr)(rep)), M, mc_iters0, mc_iters, seed + rep);
  }
}
}
