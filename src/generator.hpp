#include "utils.hpp"

class Generator
{
public:
  Generator(potts_model, int, int, std::string);
  void run(int, int);
  void writeAASequences(std::string);

private:
  int N;         // number of positions
  int Q;         // number of amino acids
  int runs;      // number of independent sampling runs
  int run_count; // number of sequences sampled from independent runs
  int M;         // number of sequences (runs * run_count)

  long int random_seed;
  int t_wait;
  int delta_t;
  double temperature;

  arma::Cube<int> samples;
  potts_model model;

  void loadParameters(std::string);
  void initializeParameters(void);
  void setParameter(std::string, std::string);

  char convertAA(int);
};
