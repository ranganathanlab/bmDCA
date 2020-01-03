#ifndef MSA_STATS_HPP
#define MSA_STATS_HPP

#include "msa.hpp"

#include <armadillo>

class MSAStats
{
public:
  MSAStats(MSA);
  // arma::Mat<double> GetRelEntropyGradient();
  // arma::Mat<double> GetFrequency1p();
  double getEffectiveM();
  double getN();
  double getM();
  double getQ();
  void writeRelEntropyGradient(std::string);
  void writeFrequency1p(std::string);
  void writeFrequency2p(std::string);
  // void writeFrequency3p(std::string);

  arma::Mat<double> frequency_1p;
  arma::field<arma::Mat<double>> frequency_2p;
  // arma::field<arma::Cube<double>> frequency_3p;
  arma::Mat<double> rel_entropy_grad_1p;

private:
  double pseudocount;
  int M;              // number of sequences
  int N;              // number of positions
  int Q;              // amino acid alphabet size
  double M_effective; // effect number of sequences

  arma::Col<double> aa_background_frequencies;
  // arma::Mat<double> frequency_1p;
  // arma::field<arma::Mat<double>> frequency_2p;
  // arma::Mat<double> rel_entropy_grad_1p;
};

#endif
