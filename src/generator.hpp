#include "mcmc_stats.hpp"
#include "graph.hpp"
#include "utils.hpp"

class Generator
{
  public:
    Generator(potts_model, int, int);
    Generator(potts_model, int, int, std::string);
    void loadParameters(std::string);
    void initializeParameters(void);
    void sample(arma::Cube<int>*);

  private:
    int N;
    int Q;
    long int seed;
    int t_wait_0;
    int t_wait;
    int delta_t;
    int delta_t_0;
    double temperature;
    bool use_independent_samples = false;

    Graph graph;
    potts_model model;
};
