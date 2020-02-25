#include <armadillo>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "generator.hpp"

int
main(int argc, char* argv[])
{
  int num_sequences = 1;

  std::string parameters_file, h_file, J_file;
  std::string dest_dir = ".";
  std::string config_file;
  std::string output_file = "mcmc_sequences.fasta";

  bool dest_dir_given = false;
  bool compat_mode = true;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "i:h:j:d:n:c:o:")) != -1) {
    switch (c) {
      case 'i':
        parameters_file = optarg;
        break;
      case 'h':
        h_file = optarg;
        compat_mode = false;
        break;
      case 'j':
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
        num_sequences = std::stoi(optarg);
        break;
      case 'c':
        config_file = optarg;
        break;
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

  Generator generator = Generator(params, N, Q, config_file);

  if (dest_dir_given == true) {
    chdir(dest_dir.c_str());
  }

  generator.run(num_sequences, 1);
  generator.writeAASequences(output_file);

  return 0;
};
