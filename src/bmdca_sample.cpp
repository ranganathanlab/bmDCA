#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <armadillo>

#include "generator.hpp"

int
main(int argc, char* argv[])
{
  int num_sequences = 1;
  int reps = 1;
  // double temperature = 1.0;

  std::string parameters_file, h_file, J_file;
  std::string dest_dir = ".";
  std::string config_file;
  std::string output_file = "mcmc_sequences.txt";

  bool dest_dir_given = false;
  bool compat_mode = true;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "i:j:d:n:c:o:")) != -1) {
    switch (c) {
      case 'i':
        parameters_file = optarg;
        // input_file_given = true;
        break;
      case 'j':
        h_file = parameters_file;
        J_file = optarg;
        compat_mode = false;
        break;
      case 'd':
        dest_dir = optarg;
        {
          struct stat st = { 0 };
          if (stat(dest_dir.c_str(), &st) == -1) {
            mkdir(dest_dir.c_str(), 0700);
          }
        }
        dest_dir_given = true;
        break;
      case 'o':
        output_file = optarg;
        break;
      case 'n':
        // num_sequences = std::stoi(optarg);
        reps = std::stoi(optarg);
        break;
      case 'c':
        config_file = optarg;
        // config_file_given = true;
        break;
      // case 't':
      //   temperature = std::stod(optarg);
      //   break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  // Load Potts model
  potts_model params;
  if (compat_mode) {
    params = loadPottsModelCompat(parameters_file);
  } else {
    params = loadPottsModel(h_file, J_file);
  }

  int N = params.h.n_cols;
  int Q = params.h.n_rows;

  if (dest_dir_given == true) {
    chdir(dest_dir.c_str());
  }

  Generator generator = Generator(params, N, Q, config_file);
  arma::Cube<int> samples =
    arma::Cube<int>(num_sequences, N, reps, arma::fill::zeros);
  generator.sample(&(samples));

  MCMCStats mcmc_stats = MCMCStats(&samples, &params);
  mcmc_stats.writeSamples(output_file);

  return 0;
};
