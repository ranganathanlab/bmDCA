#include "mcmc_stats.hpp"
#include "graph.hpp"
#include "utils.hpp"

class Generator
{
  public:
    Generator(potts_model, int, int);
    Generator(potts_model, int, int, std::string);
    void sample(arma::Cube<int>*);

  private:
    int N; // number of positions
    int Q; // number of amino acids
    // int M; // number of sequences

    long int random_seed;
    double adapt_up;
    double adapt_down;
    bool check_ergo;
    int t_wait;
    int delta_t;
    double temperature;
    bool use_indep_samples;
    bool output_numerical;

    Graph graph;
    potts_model model;

    void loadParameters(std::string);
    void initializeParameters(void);
    void setParameter(std::string, std::string);
};
