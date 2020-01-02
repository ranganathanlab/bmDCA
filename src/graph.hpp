#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <armadillo>
#include <iostream>
#include <string>

#include "mvector.hpp"
#include "utils.h"

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

  // void randomize(double beta = 1.);

  // void randomize_gauss(double beta = 1.);

  void read(std::istream& is);

  std::ostream& print_distribution(std::ostream& os);

  std::ostream& sample_distribution(std::ostream& os, size_t m);

  void sample_mcmc(arma::Mat<int>* ptr,
                   size_t m,
                   size_t mc_iters0,
                   size_t mc_iters,
                   long int seed);

  std::ostream& sample_mcmc(std::ostream& os,
                            size_t m,
                            size_t mc_iters0,
                            size_t mc_iters,
                            std::string const& out_energies_name,
                            long int seed);

  std::ostream& initialize_mcmc(std::ostream& os,
                                size_t m,
                                size_t mc_iters0,
                                size_t mc_iters,
                                std::string const& out_energies_name,
                                int* initial_conf,
                                double* tot_de_record,
                                double* tot_de_record2);

  std::ostream& print_parameters(std::ostream& os);

  void print_parameters(FILE* of);
};

#endif // GRAPH_HPP
