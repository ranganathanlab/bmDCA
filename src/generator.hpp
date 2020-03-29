#include "mcmc.hpp"
#include "mcmc_stats.hpp"
#include "utils.hpp"

class Generator
{
public:
  Generator(potts_model, int, int, std::string);
  ~Generator(void);
  void run(int, int, std::string);
  void writeAASequences(std::string);
  void writeNumericalSequences(std::string);

private:
  int N;         // number of positions
  int Q;         // number of amino acids
  int M;         // number of independent sampling runs
  int count_max; // number of sequences sampled from independent runs

  int resample_max;
  long int random_seed;
  int t_wait_0;
  int delta_t_0;
  bool check_ergo;
  double adapt_up_time;
  double adapt_down_time;
  double temperature;

  arma::Cube<int> samples;
  potts_model model;

  MCMC* mcmc;
  MCMCStats* mcmc_stats;

  void loadParameters(std::string);
  void initializeParameters(void);
  void checkParameters(void);
  void setParameter(std::string, std::string);

  char convertAA(int);
};
