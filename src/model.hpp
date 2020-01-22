#ifndef MODEL_HPP
#define MODEL_HPP

#include "msa_stats.hpp"
#include "utils.hpp"

class Model
{
public:
  potts_model params;
  potts_model learning_rates;
  potts_model gradient;
  int N;
  int Q;

  Model(MSAStats, double, double);

  void writeParams(std::string);
  void writeLearningRates(std::string);
  void writeGradient(std::string);
};

#endif
