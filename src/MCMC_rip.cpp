#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "graph.hpp"

using namespace std;

void
usage(char* command, ostream& os)
{

  os << "Usage: " << command << " [options]" << endl;
  os << endl;
  os << "Montecarlo." << endl;

  os << "  options:" << endl;
  os << "    -n <n>    number of nodes" << endl;
  os << "    -q <q>    number of states" << endl;
  os << "    -m <m>    number of samples" << endl;
  os << "    -s <s>    random seed" << endl;
  os << "    -t <t>    montecarlo iterations (default=100000)" << endl;
  os << "    -T <t>    montecarlo initial iterations (default=1000000)" << endl;
  os << "    -i        initial configurations provided on file initial_conf.txt"
     << endl;
  os << "    -u <str>  extra suffix for output files" << endl;
  os << "    -r <r>    number of restarts" << endl;
  os << "    -?        this help" << endl;
}

int
main(int argc, char** argv)
{
  int c;
  size_t rip;
  size_t n = 4;
  size_t q = 3;

  size_t m = 10;
  size_t mtot = 10;

  long int seed = 1;

  size_t mc_iters0 = 1000000;
  size_t mc_iters = 100000;
  size_t mc_rip = 1;

  bool gauge0 = false;
  bool init = false;

  string suffix("");

  while ((c = getopt(argc, argv, "n:q:m:s:b:t:T:0iu:r:?")) != -1) {
    switch (c) {
      case 'n':
        n = atoi(optarg);
        break;
      case 'q':
        q = atoi(optarg);
        break;
      case 'm':
        m = atoi(optarg);
        break;
      case 's':
        seed = atol(optarg);
        break;
      case 'T':
        mc_iters0 = atoi(optarg);
        break;
      case 't':
        mc_iters = atoi(optarg);
        break;
      case '0':
        gauge0 = true;
        break;
      case 'i':
        init = true;
        break;
      case 'u':
        suffix = string(optarg);
        break;
      case 'r':
        mc_rip = atoi(optarg);
        break;
      case '?':
        usage(argv[0], cout);
        exit(EXIT_SUCCESS);
        break;
      default:
        cerr << "Use -? for help." << endl;
        exit(EXIT_FAILURE);
    }
  }
  srand48(seed);

  mtot = m * mc_rip;

  string out_samples_name;

  stringstream ss;
  ss << "out_samples_montecarlo_n" << n << "_q" << q << "_s" << seed << "_m"
     << m << "_T" << mc_iters0 << "_t" << mc_iters << suffix << ".txt";
  ss >> out_samples_name;

  ofstream out_samples(out_samples_name.c_str());

  string out_energies_name;

  stringstream se;
  se << "out_energies_montecarlo_n" << n << "_q" << q << "_s" << seed << "_m"
     << m << "_T" << mc_iters0 << "_t" << mc_iters << suffix << ".txt";
  se >> out_energies_name;

  out_samples << mtot << " " << n << " " << q << endl;

  if (not init) {
    Graph graph(n, q);
    graph.read(cin);
    for (rip = 0; rip < mc_rip; rip++) {
      graph.sample_mcmc(
        out_samples, m, mc_iters0, mc_iters, out_energies_name, seed + rip);
    }
  } else {
    int* initial_conf;
    initial_conf = (int*)malloc(sizeof(int) * n);
    FILE* fp;
    int i;
    fp = fopen("initial_conf.txt", "r");
    int a1, a2, a3;
    fscanf(fp, "%d %d %d\n", &a1, &a2, &a3);
    double* tot_de_record;
    tot_de_record = (double*)malloc(sizeof(double) * mc_iters);
    double* tot_de_record2;
    tot_de_record2 = (double*)malloc(sizeof(double) * mc_iters);
    Graph graph(n, q);
    graph.read(cin);
    for (rip = 0; rip < mc_rip; rip++) {
      for (i = 0; i < n; i++) {
        fscanf(fp, "%d", &initial_conf[i]);
      }
      fscanf(fp, "\n");
      graph.initialize_mcmc(out_samples,
                            m,
                            mc_iters0,
                            mc_iters,
                            out_energies_name,
                            initial_conf,
                            tot_de_record,
                            tot_de_record2);
    }
    fclose(fp);
    fp = fopen("tot_de_record.txt", "w");
    for (i = 0; i < mc_iters; i++) {
      fprintf(fp,
              "%lf %lf\n",
              tot_de_record[i] / mc_rip,
              (tot_de_record2[i] / mc_rip));
    }
    fclose(fp);
    free(initial_conf);
    free(tot_de_record);
    free(tot_de_record2);
  }

  out_samples.close();

  return EXIT_SUCCESS;
}
