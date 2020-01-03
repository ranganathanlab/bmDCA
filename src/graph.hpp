#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <armadillo>
#include <iostream>
#include <string>

#include "mvector.hpp"
#include "utils.hpp"

class Graph
{
public:
  Graph(size_t n, size_t q)
    : n(n)
    , q(q)
    , J(xstd::mshape<4>(n, n, q, q))
    , h(xstd::mshape<2>(n, q)){};

  void load(potts_model);

  size_t n, q;
  xstd::mvector<4, double> J;
  xstd::mvector<2, double> h;

  std::ostream& print_distribution(std::ostream& os);

  std::ostream& print_parameters(std::ostream& os);

  // std::ostream& sample_distribution(std::ostream& os, size_t m);

  void sample_mcmc(arma::Mat<int>* ptr,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   long int seed);

  void sample_mcmc_init(arma::Mat<int>* ptr,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   arma::Col<int>* init_ptr);

  void print_parameters(FILE* of);
};

#endif
