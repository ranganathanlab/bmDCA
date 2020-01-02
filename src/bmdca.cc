#include <armadillo>

#include "msa.h"
#include "msa_stats.h"
#include "run.h"

int
main(int argc, char* argv[])
{
  std::string input_file;
  std::string dest_dir;
  bool reweight = false;

  // Read command-line parameters
  char c;
  while ((c = getopt(argc, argv, "i:d:r")) != -1) {
    switch (c) {
      case 'i':
        input_file = optarg;
        break;
      case 'd':
        dest_dir = optarg;
        break;
      case 'r':
        reweight = true;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  MSA msa = MSA(input_file, reweight);
  MSAStats msa_stats = MSAStats(msa);
  Sim sim = Sim(msa_stats);
  // sim.load(params);
  sim.run();

  return 0;
}
