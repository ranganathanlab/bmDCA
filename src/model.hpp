#ifndef MODEL_HPP
#define MODEL_HPP

#include "msa_stats.hpp"
#include "utils.hpp"

class Model
{
public:
  potts_model params;
  potts_model params_prev;
  potts_model learning_rates;
  potts_model gradient;
  potts_model gradient_prev;
  int N;
  int Q;

  Model(MSAStats, double, double);
  Model(std::string, std::string, std::string, std::string, std::string);
  Model(std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string,
        std::string);

  void writeParams(std::string, std::string);
  void writeParamsPrevious(std::string, std::string);
  void writeLearningRates(std::string, std::string);
  void writeGradient(std::string, std::string);
  void writeGradientPrevious(std::string, std::string);

  void writeParamsCompat(std::string);
  void writeParamsPreviousCompat(std::string);
  void writeLearningRatesCompat(std::string);
  void writeGradientCompat(std::string);
  void writeGradientPreviousCompat(std::string);
};

#endif
