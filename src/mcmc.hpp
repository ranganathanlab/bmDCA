#ifndef MCMC_HPP
#define MCMC_HPP

#include <string>
#include <unistd.h>

#include "graph.hpp"
#include "utils.hpp"

class MCMC
{

public:
  MCMC(size_t N, size_t Q);
  MCMC(potts_model, size_t N, size_t Q);
  void load(potts_model);
  void run(int, int);
  void
  sample(arma::Cube<int>*, int, int, int, int, int, long int, double);
  void sample_init(arma::Cube<int>*,
                   int,
                   int,
                   int,
                   int,
                   int,
                   arma::Col<int>*,
                   double);

private:
  size_t n; // number of positions
  size_t q; // number of amino acids (inc. gaps)

  Graph graph;
};

#endif
