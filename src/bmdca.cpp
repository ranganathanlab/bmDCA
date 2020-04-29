#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "msa.hpp"
#include "msa_stats.hpp"
#include "run.hpp"

void
print_usage(void)
{
  std::cout << "bmdca usage:" << std::endl;
  std::cout << "(e.g. bmdca -i <input MSA> -r -d <directory> -c <config file>)"
            << std::endl;
  std::cout << "  -i: input MSA (FASTA format)" << std::endl;
  std::cout << "  -d: destination directory" << std::endl;
  std::cout << "  -r: re-weighting flag" << std::endl;
  std::cout << "  -n: numerical MSA" << std::endl;
  std::cout << "  -w: sequence weights" << std::endl;
  std::cout << "  -c: config file" << std::endl;
  std::cout << "  -t: sequence similarity threshold for reweighting" << std::endl;
  std::cout << "  -h: print usage (i.e. this message)" << std::endl;
  std::cout << "  -f: force a restart of the inference loop" << std::endl;
};

int
main(int argc, char* argv[])
{
  std::string input_file;
  std::string config_file;
  std::string weight_file;
  std::string dest_dir = ".";

  bool reweight = false; // if true, weight sequences by similarity
  bool numeric_msa_given = false;
  bool input_file_given = false;
  bool force_restart = false;
  bool weight_given = false;
  double threshold = 0.8;

  // Read command-line parameters.
  char c;
  while ((c = getopt(argc, argv, "c:d:fhi:n:rt:w:")) != -1) {
    switch (c) {
      case 'c':
        config_file = optarg;
        break;
      case 'd':
        dest_dir = optarg;
        {
          struct stat st = { 0 };
          if (stat(dest_dir.c_str(), &st) == -1) {
#if __unix__
            mkdir(dest_dir.c_str(), 0700);
#elif _WIN32
            mkdir(dest_dir.c_str());
#endif
          }
        }
        break;
      case 'f':
        force_restart = true;
        break;
      case 'h':
        print_usage();
        return 0;
        break;
      case 'i':
        input_file = optarg;
        input_file_given = true;
        break;
      case 'n':
        input_file = optarg;
        input_file_given = true;
        numeric_msa_given = true;
        break;
      case 'r':
        reweight = true;
        break;
      case 't':
        threshold = std::stod(optarg);
        break;
      case 'w':
        weight_file = optarg;
        weight_given = true;
        break;
      case '?':
        std::cerr << "ERROR: Incorrect command line usage." << std::endl;
        print_usage();
        std::exit(EXIT_FAILURE);
    }
  }

  // Exit if no input file was provided.
  if (!input_file_given) {
    std::cerr << "ERROR: Missing MSA input file." << std::endl;
    print_usage();
    std::exit(EXIT_FAILURE);
  }

  // If both the numeric matrix and sequence weights are given, don't bother
  // converting the FASTA file.
  if (numeric_msa_given && weight_given) {
    // Parse the multiple sequence alignment. Reweight sequences if desired.
    MSA msa = MSA(input_file, weight_file, numeric_msa_given);
    msa.writeSequenceWeights(dest_dir + "/sequence_weights.txt");
    msa.writeMatrix(dest_dir + "/msa_numerical.txt");

    // Compute the statistics of the MSA.
    MSAStats msa_stats = MSAStats(msa);
    msa_stats.writeFrequency1p(dest_dir + "/stat_align_1p.bin");
    msa_stats.writeFrequency2p(dest_dir + "/stat_align_2p.bin");
    msa_stats.writeRelEntropyGradient(dest_dir + "/rel_ent_grad_align_1p.txt");

    // Initialize the MCMC using the statistics of the MSA.
    Sim sim = Sim(msa_stats, config_file, dest_dir, force_restart);
    sim.writeParameters("bmdca_params.conf");
    sim.run();
  } else {
    // Parse the multiple sequence alignment. Reweight sequences if desired.
    MSA msa = MSA(input_file, reweight, numeric_msa_given, threshold);
    msa.writeSequenceWeights(dest_dir + "/sequence_weights.txt");
    msa.writeMatrix(dest_dir + "/msa_numerical.txt");

    // Compute the statistics of the MSA.
    MSAStats msa_stats = MSAStats(msa);
    msa_stats.writeFrequency1p(dest_dir + "/stat_align_1p.bin");
    msa_stats.writeFrequency2p(dest_dir + "/stat_align_2p.bin");
    msa_stats.writeRelEntropyGradient(dest_dir + "/rel_ent_grad_align_1p.txt");

    // Initialize the MCMC using the statistics of the MSA.
    Sim sim = Sim(msa_stats, config_file, dest_dir, force_restart);
    sim.writeParameters("bmdca_params.conf");
    sim.run();
  }

  return 0;
};
