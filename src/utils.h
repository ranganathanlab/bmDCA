#ifndef UTILS_H
#define UTILS_H

#include <armadillo>
#include <string>

class SeqRecord
{
private:
  const std::string header;
  const std::string sequence;

public:
  SeqRecord(std::string, std::string);
  // virtual ~SeqRecord(void);
  void PrintRecord();
  std::string GetRecord();
  std::string GetHeader();
  std::string GetSequence();
};

typedef struct
{
  arma::field<arma::Mat<double>> J;
  arma::Mat<double> h;
  // int N; // number of positions
  // int Q; // number of amino acids (21, including gaps)
} potts_model;

void WritePottsModelCompat(potts_model, std::string);

void WriteMCMCSamplesCompat(arma::field<arma::Mat<int>>, std::string);

void WriteMCMCEnergiesCompat(arma::Col<double>, std::string);

double
maximum(double, double);

int
Theta(double x);
int
Delta(double x);
double
Max(double a, double b);
double
Min(double a, double b);

#endif
