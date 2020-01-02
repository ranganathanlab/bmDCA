#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "msa.h"
#include "msa_stats.h"

int
main(int argc, char* argv[])
{

  std::string infile;
  std::string dest;
  bool reweight = false;

  char c;
  while ((c = getopt(argc, argv, "i:d:r")) != -1) {
    switch (c) {
      case 'i':
        infile = optarg;
        break;
      case 'd':
        dest = optarg;
        break;
      case 'r':
        reweight = true;
        break;
      case '?':
        std::cerr << "what the fuck?" << std::endl;
    }
  }

  MSA msa = MSA(infile, reweight);
  msa.WriteSequenceWeightsCompat(dest + "/weights.txt");
  msa.WriteMatrixCompat(dest + "/msa_numerical.txt");
  MSAStats msa_stats = MSAStats(msa);
  msa_stats.WriteFrequency1pCompat(dest + "/stat_align_1p.txt");
  msa_stats.WriteFrequency2pCompat(dest + "/stat_align_2p.txt");
  // msa_stats.WriteFrequency3pCompat(dest + "/stat_align_3p.txt");
  msa_stats.WriteRelEntropyGradientCompat(dest +
                                          "/relentropy_grad_align_1p.txt");
  return 0;
}
