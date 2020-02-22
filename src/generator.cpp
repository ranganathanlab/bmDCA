#include "generator.hpp"

Generator::Generator(potts_model params,
                                     int n,
                                     int q,
                                     std::string config_file)
  : N(n), Q(q), graph(n, q)
{
  if (config_file.length() == 0) {
    initializeParameters();
  } else {
    loadParameters(config_file);
  }
  graph.load(params);
};

void
Generator::initializeParameters(void) {
  seed = 1;
  t_wait = 100000;
  delta_t = 1000;
  temperature = 1.0;
  use_independent_samples = true;
};

void
Generator::loadParameters(std::string config_file) {
}

void
Generator::sample(arma::Cube<int>* ptr) {
  int M = ptr->n_rows;
  int N = ptr->n_cols;
  int reps = ptr->n_slices;

  if (use_independent_samples) {
#pragma omp parallel
    {
#pragma omp for
      for (int rep = 0; rep < reps; rep++) {
        graph.sample_mcmc((arma::Mat<int>*)&(*ptr).slice(rep),
                          M,
                          t_wait,
                          delta_t,
                          seed + rep,
                          temperature);
      }
    }
  }
};

