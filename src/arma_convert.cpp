#include <armadillo>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "utils.hpp"

void
print_usage(void)
{
  std::cout << "arma2ascii usage:" << std::endl;
  std::cout << "(e.g. bmdca_sample -p <params h> -P <params J>" << std::endl;
  std::cout << " -OR- bmdca_sample -p <params>" << std::endl;
  std::cout << " -OR- bmdca_sample -s <stats file>)" << std::endl;
  std::cout << "  -p: parameters (txt) _or_ fields h (bin)" << std::endl;
  std::cout << "  -P: couplings J (bin), *required* if fields h given"
            << std::endl;
  std::cout << "  -s: sequence sample statistics file" << std::endl;
  std::cout << "  -h: print usage (i.e. this message)" << std::endl;
};

int
main(int argc, char* argv[])
{
  std::string stat_file;
  std::string param_h_file;
  std::string param_J_file;

  bool valid_input = false;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "P:p:s:h")) != -1) {
    switch (c) {
      case 'p':
        param_h_file = optarg;
        valid_input = true;
        break;
      case 'P':
        param_J_file = optarg;
        break;
      case 's':
        stat_file = optarg;
        valid_input = true;
        break;
      case 'h':
        print_usage();
        return 0;
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        print_usage();
        std::exit(EXIT_FAILURE);
    }
  }

  if (!valid_input) {
    std::cerr << "ERROR: Input file not given." << std::endl;
    print_usage();
    std::exit(EXIT_FAILURE);
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
