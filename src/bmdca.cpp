// #include <armadillo>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "msa.hpp"
#include "msa_stats.hpp"
#include "run.hpp"

int
main(int argc, char* argv[])
{
  std::string input_file;
  std::string config_file;
  std::string dest_dir;
  bool reweight = false;
  bool dest_dir_given = false;

  // Read command-line parameters
  char c;
  while ((c = getopt(argc, argv, "i:d:c:r")) != -1) {
    switch (c) {
      case 'i':
        input_file = optarg;
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
      case 'c':
        config_file = optarg;
        break;
      case 'r':
        reweight = true;
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
    }
  }

  // Parse the multiple sequence alignment. Reweight sequences if desired.
  MSA msa = MSA(input_file, reweight);
  msa.writeSequenceWeights(dest_dir + "/sequence_weights.txt");
  msa.writeMatrix(dest_dir + "/msa_numerical.txt");

  // Compute the statistics of the MSA
  MSAStats msa_stats = MSAStats(msa);
  msa_stats.writeFrequency1p(dest_dir + "/stat_align_1p.txt");
  msa_stats.writeFrequency2p(dest_dir + "/stat_align_2p.txt");
  msa_stats.writeRelEntropyGradient(dest_dir + "/rel_ent_grad_align_1p.txt");

  // Initialize the MCMC using the statistics of the MSA
  Sim sim = Sim(msa_stats);
  // sim.load(params);

  if (dest_dir_given == true) {
    chdir(dest_dir.c_str());
  }

  sim.run();

  return 0;
}
