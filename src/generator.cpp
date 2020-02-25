#include "generator.hpp"

#include "mcmc.hpp"
// #include "mcmc_stats.hpp"

Generator::Generator(potts_model params, int n, int q, std::string config_file)
  : N(n)
  , Q(q)
  , model(params)
{
  if (config_file.length() == 0) {
    initializeParameters();
  } else {
    loadParameters(config_file);
  }
};

void
Generator::initializeParameters(void)
{
  random_seed = 1;
  t_wait = 10000;
  delta_t = 100; // placeholder
  temperature = 1.0;
};

void
Generator::loadParameters(std::string file_name)
{
  std::ifstream file(file_name);
  bool reading_bmdca_section = false;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '#' || line.empty()) {
        continue;
      } else if (line[0] == '[') {
        if (line == "[sampling]") {
          reading_bmdca_section = true;
          continue;
        } else {
          reading_bmdca_section = false;
          continue;
        }
      }
      if (reading_bmdca_section) {
        auto delim_pos = line.find("=");
        auto key = line.substr(0, delim_pos);
        auto value = line.substr(delim_pos + 1);
        setParameter(key, value);
      }
    }
  } else {
    std::cerr << "ERROR: " << file_name << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }
};

void
Generator::setParameter(std::string key, std::string value)
{
  if (key == "random_seed") {
    random_seed = std::stoi(value);
  } else if (key == "t_wait") {
    t_wait = std::stoi(value);
  } else if (key == "temperature") {
    temperature = std::stod(value);
  }
};

void
Generator::writeAASequences(std::string output_file)
{
  int M = samples.n_rows;
  int N = samples.n_cols;
  int reps = samples.n_slices;

  std::ofstream output_stream(output_file);

  char aa;
  for (int rep = 0; rep < reps; rep++) {
    for (int m = 0; m < M; m++) {
      output_stream << ">sample" << m * rep + rep << std::endl;
      for (int n = 0; n < N; n++) {
        aa = convertAA(samples.at(m, n, rep));
        if (aa != '\0') {
          output_stream << aa;
        }
      }
      output_stream << std::endl << std::endl;
    }
  }
};

char
Generator::convertAA(int n)
{
  char aa = '\0';
  switch (n) {
    case 0:
      break;
    case 1:
      aa = 'A';
      break;
    case 2:
      aa = 'C';
      break;
    case 3:
      aa = 'D';
      break;
    case 4:
      aa = 'E';
      break;
    case 5:
      aa = 'F';
      break;
    case 6:
      aa = 'G';
      break;
    case 7:
      aa = 'H';
      break;
    case 8:
      aa = 'I';
      break;
    case 9:
      aa = 'K';
      break;
    case 10:
      aa = 'L';
      break;
    case 11:
      aa = 'M';
      break;
    case 12:
      aa = 'N';
      break;
    case 13:
      aa = 'P';
      break;
    case 14:
      aa = 'Q';
      break;
    case 15:
      aa = 'R';
      break;
    case 16:
      aa = 'S';
      break;
    case 17:
      aa = 'T';
      break;
    case 18:
      aa = 'V';
      break;
    case 19:
      aa = 'W';
      break;
    case 20:
      aa = 'Y';
      break;
  }
  return aa;
};

void
Generator::run(int n_indep_runs, int n_per_run)
{
  runs = n_indep_runs;
  run_count = n_per_run;
  M = runs * run_count;

  samples = arma::Cube<int>(run_count, N, runs, arma::fill::zeros);

  MCMC mcmc = MCMC(model, N, Q);
  mcmc.sample(
    &samples, runs, run_count, N, t_wait, delta_t, random_seed, temperature);

  // MCMCStats mcmc_stats = MCMCStats(&samples, &model);
  // mcmc_stats.writeSamples("MC_samples.txt");
  // mcmc_stats.writeSampleEnergies("MC_energies.txt");
};
