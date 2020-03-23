#include "utils.hpp"

#include <cstdio>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/types.h>

#ifndef AA_ALPHABET_SIZE
#define AA_ALPHABET_SIZE 21
#endif

SeqRecord::SeqRecord(std::string h, std::string s)
  : header(h)
  , sequence(s){};

void
SeqRecord::print(void)
{
  std::cout << ">" << header << std::endl;
  std::cout << sequence << std::endl;
};

std::string
SeqRecord::getRecord(void)
{
  std::string record_string = ">" + header + "\n" + sequence + "\n";
  return record_string;
};

std::string
SeqRecord::getHeader(void)
{
  return header;
};

std::string
SeqRecord::getSequence(void)
{
  return sequence;
};

potts_model
loadPottsModel(std::string h_file, std::string J_file)
{
  potts_model params;
  params.h.load(h_file);
  params.J.load(J_file);
  return params;
};

potts_model
loadPottsModelHDF(std::string parameters_file)
{
  arma::Mat<int> N_mat = arma::Mat<int>(1, 1);
  arma::Mat<int> Q_mat = arma::Mat<int>(1, 1);

  N_mat.load(arma::hdf5_name(parameters_file, "N"));
  Q_mat.load(arma::hdf5_name(parameters_file, "Q"));

  int N = (int)as_scalar(N_mat);
  int Q = (int)as_scalar(Q_mat);

  potts_model params;
  params.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  params.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J.at(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

  params.h.load(arma::hdf5_name(parameters_file, "h"));

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J.at(i, j).load(arma::hdf5_name(
        parameters_file, "J/" + std::to_string(i) + "_" + std::to_string(j)));
    }
  }
  return params;
};

potts_model
loadPottsModelCompat(std::string parameters_file)
{
  std::ifstream input_stream(parameters_file);

  if (!input_stream) {
    std::cerr << "ERROR: couldn't open '" << parameters_file << "' for reading."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int count = 0;
  std::string line;
  while (std::getline(input_stream, line))
    count++;

  int N;
  int Q = AA_ALPHABET_SIZE;
  N =
    (int)(sqrt(2 * count + ((double)Q - 2) * ((double)Q - 2) / 4) / (double)Q +
          ((double)Q - 2) / (2 * (double)Q));

  input_stream.clear();
  input_stream.seekg(0);

  potts_model params;
  params.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  params.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J.at(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

  // Read parameters
  int n1, n2, aa1, aa2;
  double value;
  std::string tmp = "";
  for (int count = 0; count < (int)N * (N - 1) / 2 * Q * Q; count++) {
    input_stream >> tmp;
    input_stream >> n1 >> n2 >> aa1 >> aa2;
    input_stream >> value;
    params.J.at(n1, n2).at(aa1, aa2) = value;
  }

  for (int count = 0; count < N * Q; count++) {
    input_stream >> tmp;
    input_stream >> n1 >> aa1;
    input_stream >> value;
    params.h.at(aa1, n1) = value;
  }
  return params;
};

int
Theta(double x)
{
  if (x > 0)
    return 1;
  return 0;
};

int
Delta(double x)
{
  if (x == 0)
    return 1;
  return 0;
};

double
Max(double a, double b)
{
  if (a > b)
    return a;
  return b;
};

double
Min(double a, double b)
{
  if (a < b)
    return a;
  return b;
};

int
deleteFile(std::string filename) {
  std::fstream fs;
  fs.open(filename);
  if (!fs.fail()) {
    if (std::remove(filename.c_str()) != 0) {
      fs.close();
      return 1;
    }
  }
  fs.close();
  return 0;
};

bool
checkFileExists(std::string filename) {
  std::fstream fs;
  fs.open(filename);
  if (fs.fail()) {
    return false;
  } else {
    return true;
  }
};

void
deleteAllFiles(std::string directory) {
  DIR *dp;
  struct dirent *dirp;

  dp = opendir(".");
  std::vector<int> steps;
  while ((dirp = readdir(dp)) != NULL) {
    std::string fname = dirp->d_name;
    if (fname == ".") continue;
    if (fname == "..") continue;
    if (std::remove(fname.c_str()) != 0) {
      std::cerr << "ERROR: deletion of '" << fname << "' failed." << std::endl;
    }
  }
  closedir(dp);
  return;
};
