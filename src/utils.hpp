#ifndef UTILS_HPP
#define UTILS_HPP

#include <armadillo>
#include <string>

class SeqRecord
{
private:
  const std::string header;
  const std::string sequence;

public:
  SeqRecord(std::string, std::string);
  void print();
  std::string getRecord();
  std::string getHeader();
  std::string getSequence();
};

typedef struct
{
  arma::field<arma::Mat<double>> J;
  arma::Mat<double> h;
} potts_model;

int
Theta(double);

int
Delta(double);

double
Max(double, double);

double
Min(double, double);

#endif
