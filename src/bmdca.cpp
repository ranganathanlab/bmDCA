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
  std::string numeric_file;
  std::string weight_file;
  std::string dest_dir = ".";

  bool reweight = false;
  bool dest_dir_given = false;
  bool numeric_msa_given = false;
  bool input_file_given = true;
  bool weight_given = false;
  double threshold = 0.8;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "i:d:c:rn:w:t:")) != -1) {
    switch (c) {
      case 'i':
        input_file = optarg;
        input_file_given = true;
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
      case 'n':
        numeric_file = optarg;
        numeric_msa_given = true;
        break;
      case 'w':
        weight_file = optarg;
        weight_given = true;
        break;
      case 't':
        threshold = std::stod(optarg);
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        std::exit(EXIT_FAILURE);
    }
  }

  // If both the numeric matrix and sequence weights are given, don't bother
  // converting the FASTA file.
  if (numeric_msa_given && weight_given) {
    // Parse the multiple sequence alignment. Reweight sequences if desired.
    MSA msa = MSA(numeric_file, weight_file, numeric_msa_given);
    msa.writeSequenceWeights(dest_dir + "/sequence_weights.txt");
    msa.writeMatrix(dest_dir + "/msa_numerical.txt");

    // Compute the statistics of the MSA.
    MSAStats msa_stats = MSAStats(msa);
    msa_stats.writeFrequency1p(dest_dir + "/stat_align_1p.txt");
    msa_stats.writeFrequency2p(dest_dir + "/stat_align_2p.txt");
    msa_stats.writeRelEntropyGradient(dest_dir + "/rel_ent_grad_align_1p.txt");

    // Initialize the MCMC using the statistics of the MSA.
    Sim sim = Sim(msa_stats, config_file);

    if (dest_dir_given == true) {
      chdir(dest_dir.c_str());
    }

    sim.writeParameters("bmdca_params.conf");
    sim.run();
  } else if (numeric_msa_given) {
    MSA msa = MSA(numeric_file, reweight, numeric_msa_given, threshold);
    msa.writeSequenceWeights(dest_dir + "/sequence_weights.txt");
    msa.writeMatrix(dest_dir + "/msa_numerical.txt");

    // Compute the statistics of the MSA.
    MSAStats msa_stats = MSAStats(msa);
    msa_stats.writeFrequency1p(dest_dir + "/stat_align_1p.txt");
    msa_stats.writeFrequency2p(dest_dir + "/stat_align_2p.txt");

    // Initialize the MCMC using the statistics of the MSA.
    Sim sim = Sim(msa_stats, config_file);

    if (dest_dir_given == true) {
      chdir(dest_dir.c_str());
    }

    sim.writeParameters("bmdca_params.conf");
    sim.run();
  } else if (input_file_given) {
    // Parse the multiple sequence alignment. Reweight sequences if desired.
    MSA msa = MSA(input_file, reweight, numeric_msa_given, threshold);
    msa.writeSequenceWeights(dest_dir + "/sequence_weights.txt");
    msa.writeMatrix(dest_dir + "/msa_numerical.txt");

    // Compute the statistics of the MSA.
    MSAStats msa_stats = MSAStats(msa);
    msa_stats.writeFrequency1p(dest_dir + "/stat_align_1p.txt");
    msa_stats.writeFrequency2p(dest_dir + "/stat_align_2p.txt");
    msa_stats.writeRelEntropyGradient(dest_dir + "/rel_ent_grad_align_1p.txt");

    // Initialize the MCMC using the statistics of the MSA.
    Sim sim = Sim(msa_stats, config_file);

    if (dest_dir_given == true) {
      chdir(dest_dir.c_str());
    }

    sim.writeParameters("bmdca_params.conf");
    sim.run();
  }

  return 0;
};
