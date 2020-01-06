#ifndef MCMC_STATS_HPP
#define MCMC_STATS_HPP

#include <armadillo>

#include "utils.hpp"

class MCMCStats
{
public:
  MCMCStats(arma::field<arma::Mat<int>>, potts_model);
  void updateData(arma::field<arma::Mat<int>>, potts_model);

  void computeEnergies(void);
  void computeEnergiesStats(void);
  void computeAutocorrelation(void);
  void computeCorrelationStats(void);
  void computeSampleStats(void);
  void computeSampleStatsImportance(potts_model* cur, potts_model* prev);

  std::vector<double> getEnergiesStats(void);
  std::vector<double> getCorrelationsStats(void);

  void writeEnergyStats(std::string, std::string, std::string, std::string);
  void writeAutocorrelationStats(std::string, std::string, std::string);
  void writeFrequency1p(std::string, std::string);
  void writeFrequency2p(std::string, std::string);

  void writeSamples(std::string);
  void writeSampleEnergies(std::string);
  void writeSampleEnergiesRelaxation(std::string, int=1);

  arma::Mat<double> frequency_1p;
  arma::Mat<double> frequency_1p_sigma;
  arma::field<arma::Mat<double>> frequency_2p;
  arma::field<arma::Mat<double>> frequency_2p_sigma;

private:
  potts_model params;
  arma::field<arma::Mat<int>> samples;
  arma::Mat<double> energies;

  double energies_start_avg;
  double energies_start_sigma;
  double energies_end_avg;
  double energies_end_sigma;
  double energies_err;
  arma::Row<double> energies_relax;
  arma::Row<double> energies_relax_sigma;

  arma::Col<double> overlaps;
  arma::Col<double> overlaps_sigma;

  double overlap_inf;
  double overlap_inf_sigma;
  double overlap_auto;
  double overlap_cross;
  double overlap_check;
  double sigma_auto;
  double sigma_cross;
  double sigma_check;
  double err_cross_auto;
  double err_cross_check;
  double err_check_auto;

  int reps;
  int N;
  int Q;
  int M;
};

#endif
