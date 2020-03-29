#include <armadillo>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "utils.hpp"

int
main(int argc, char* argv[])
{
  std::string stat_file;
  std::string param_h_file;
  std::string param_J_file;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "P:p:s:")) != -1) {
    switch (c) {
      case 'p':
        param_h_file = optarg;
        break;
      case 'P':
        param_J_file = optarg;
        break;
      case 's':
        stat_file = optarg;
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  // Convert 1p and/or 2p statistics from arma binary to ascii.
  if (stat_file.size() != 0) {
    std::cout << "converting '" << stat_file << "' to text... " << std::flush;
    convertFrequencyToAscii(stat_file);
    std::cout << "done" << std::endl;
  }

  // Convert h and J parameters in binary format to one text file.
  if ((param_h_file.size() != 0) & (param_J_file.size() != 0)) {
    std::cout << "converting '" << param_h_file << "' and '" << param_J_file
              << "' to text... " << std::flush;
    convertParametersToAscii(param_h_file, param_J_file);
    std::cout << "done" << std::endl;
  } else if ((param_h_file.size() != 0) || (param_J_file.size() != 0)) {
    std::cerr << "parameters file missing. exiting..." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return 0;
}
