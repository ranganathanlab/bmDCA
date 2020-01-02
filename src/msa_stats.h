#ifndef MSA_STATS_H
#define MSA_STATS_H

#include "msa.h"

#include <armadillo>

class MSAStats
{
public:
  MSAStats(MSA);
  arma::Mat<double> GetRelEntropyGradient();
  void WriteRelEntropyGradientCompat(std::string);
  void WriteFrequency1pCompat(std::string);
  void WriteFrequency2pCompat(std::string);

private:
  double pseudocount;
  bool reweight;
  int M;
  int L;
  double M_effective;

  arma::Col<double> aa_background_frequencies;
  arma::Mat<double> frequency_1p;
  arma::field<arma::Mat<double>> frequency_2p;
  arma::Mat<double> rel_entropy_grad_1p;
};

#endif
