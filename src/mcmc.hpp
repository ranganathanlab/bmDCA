#ifndef MCMC_HPP
#define MCMC_HPP

#include <unistd.h>
#include <string>

#include "graph.hpp"
#include "utils.hpp"

class MCMC
{

public:
  MCMC(size_t N, size_t Q);
  MCMC(potts_model, size_t N, size_t Q);
  void load(potts_model);
  void run(int, int);
  void sample(arma::field<arma::Mat<int>>*, int, int, int, int, int, int);
  void sample_init(arma::field<arma::Mat<int>>*, int, int, int, int, int, arma::Col<int>*);

private:
  size_t n;
  size_t q;
  
  long int seed = 1;

  Graph graph;
};

#endif
