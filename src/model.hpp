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
  Model(std::string, std::string, std::string);
  Model(std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string);

  void writeParams(std::string, std::string);
  void writeLearningRates(std::string, std::string);
  void writeGradient(std::string, std::string);

  void writeParamsCompat(std::string);
  void writeLearningRatesCompat(std::string);
  void writeGradientCompat(std::string);
};

#endif
