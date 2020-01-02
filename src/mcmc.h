#ifndef MCMC_H
#define MCMC_H

// #include <cstdlib>
#include <unistd.h>

#include "graph.hpp"
#include "utils.h"

class MCMC
{

public:
  // MCMC(void);
  // MCMC(int N, int Q);
  MCMC(size_t N, size_t Q);
  MCMC(potts_model, size_t N, size_t Q);
  // MCMC(potts_model, size_t N, size_t Q);
  void load(potts_model);
  void run(int, int);
  void sample(arma::field<arma::Mat<int>>*, int, int, int, int, int, int);
  // void initializeMCMC(void);
  // void initializeGraph(potts_model);

private:
  size_t rip;
  size_t n;
  size_t q;

  size_t m;
  size_t mtot;

  long int seed;

  size_t mc_iters0;
  size_t mc_iters;
  size_t mc_rip;

  bool gauge0;
  // bool init;

  std::string suffix;

  Graph graph;
};

#endif
