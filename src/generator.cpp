#include "generator.hpp"

#include "mcmc.hpp"
#include "mcmc_stats.hpp"
#include "pcg_random.hpp"
#include "utils.hpp"

Generator::Generator(potts_model params, int n, int q, std::string config_file)
  : N(n)
  , Q(q)
  , model(params)
{
  initializeParameters();
  if (config_file.length() != 0) {
    loadParameters(config_file);
  }
};

Generator::~Generator(void)
{
  delete mcmc_stats;
};

void
Generator::initializeParameters(void)
{
  resample_max = 20;
  random_seed = 1;
  t_wait_0 = 100000;
  delta_t_0 = 1000;
  check_ergo = true;
  adapt_up_time = 1.5;
  adapt_down_time = 0.600;
  temperature = 1.0;
};

void
Generator::checkParameters(void)
{
  // Ensure that the set of ergodiciy checks is disabled if M=1
  if ((M == 1) && check_ergo) {
    check_ergo = false;
    std::cerr << "WARNING: disabling 'check_ergo' when M=1." << std::endl;
  }
}

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
  if (key == "resample_max") {
    resample_max = std::stoi(value);
  } else if (key == "random_seed") {
    random_seed = std::stoi(value);
  } else if (key == "t_wait_0") {
    t_wait_0 = std::stoi(value);
  } else if (key == "delta_t_0") {
    delta_t_0 = std::stoi(value);
  } else if (key == "check_ergo") {
    if (value.size() == 1) {
      check_ergo = (std::stoi(value) == 1);
    } else {
      check_ergo = (value == "true");
    }
  } else if (key == "adapt_up_time") {
    adapt_up_time = std::stod(value);
  } else if (key == "adapt_down_time") {
    adapt_down_time = std::stod(value);
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
      output_stream << ">sample" << m * rep + m << std::endl;
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

void
Generator::writeNumericalSequences(std::string output_file)
{
  int idx = output_file.find_last_of(".");
  std::string raw_file = output_file.substr(0, idx);

  mcmc_stats->writeSamples(raw_file + "_numerical.txt");
  mcmc_stats->writeSampleEnergies(raw_file + "_energies.txt");
}

char
Generator::convertAA(int n)
{
  char aa = '\0';
  switch (n) {
    case 0:
      aa = '-';
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
Generator::run(int n_indep_runs, int n_per_run, std::string output_file)
{
  std::cout << "initializing sampler... " << std::flush;

  arma::wall_clock timer;
  timer.tic();

  count_max = n_indep_runs;
  M = n_per_run;

  checkParameters();

  samples = arma::Cube<int>(M, N, count_max, arma::fill::zeros);
  mcmc = new MCMC(model, N, Q);
  mcmc->load(model);
  mcmc_stats = new MCMCStats(&samples, &(model));

  // Instantiate the PCG random number generator and unifrom random
  // distribution.
  pcg32 rng;
  rng.seed(random_seed);
  std::uniform_int_distribution<long int> dist(0, RAND_MAX - count_max);

  std::cout << timer.toc() << " sec" << std::endl << std::endl;

  int t_wait = t_wait_0;
  int delta_t = delta_t_0;
  bool flag_mc = true;
  int resample_counter = 0;
  while (flag_mc) {
    std::cout << "sampling model with mcmc... " << std::flush;
    timer.tic();
    mcmc->sample(
      &samples, count_max, M, N, t_wait, delta_t, dist(rng), temperature);
    std::cout << timer.toc() << " sec" << std::endl;

    std::cout << "updating mcmc stats with samples... " << std::flush;
    timer.tic();
    mcmc_stats->updateData(&samples, &model);
    std::cout << timer.toc() << " sec" << std::endl;

    if (check_ergo) {
      std::cout << "computing sequence energies and correlations... "
                << std::flush;
      timer.tic();
      mcmc_stats->computeEnergiesStats();
      mcmc_stats->computeCorrelations();
      std::cout << timer.toc() << " sec" << std::endl;

      std::vector<double> energy_stats = mcmc_stats->getEnergiesStats();
      std::vector<double> corr_stats = mcmc_stats->getCorrelationsStats();

      double auto_corr = corr_stats.at(2);
      double cross_corr = corr_stats.at(3);
      double check_corr = corr_stats.at(4);
      double cross_check_err = corr_stats.at(9);
      double auto_cross_err = corr_stats.at(8);

      double e_start = energy_stats.at(0);
      double e_start_sigma = energy_stats.at(1);
      double e_end = energy_stats.at(2);
      double e_end_sigma = energy_stats.at(3);
      double e_err = energy_stats.at(4);

      bool flag_deltat_up = true;
      bool flag_deltat_down = true;
      bool flag_twaiting_up = true;
      bool flag_twaiting_down = true;

      if (check_corr - cross_corr < cross_check_err) {
        flag_deltat_up = false;
      }
      if (auto_corr - cross_corr > auto_cross_err) {
        flag_deltat_down = false;
      }

      if (e_start - e_end < 2 * e_err) {
        flag_twaiting_up = false;
      }
      if (e_start - e_end > -2 * e_err) {
        flag_twaiting_down = false;
      }

      if (flag_deltat_up) {
        delta_t = (int)(round((double)delta_t * adapt_up_time));
        std::cout << "increasing wait time to " << delta_t << std::endl;
      } else if (flag_deltat_down) {
        delta_t = Max((int)(round((double)delta_t * adapt_down_time)), 1);
        std::cout << "decreasing wait time to " << delta_t << std::endl;
      }

      if (flag_twaiting_up) {
        t_wait = (int)(round((double)t_wait * adapt_up_time));
        std::cout << "increasing burn-in time to " << t_wait << std::endl;
      } else if (flag_twaiting_down) {
        t_wait = Max((int)(round((double)t_wait * adapt_down_time)), 1);
        std::cout << "decreasing burn-in time to " << t_wait << std::endl;
      }

      if (not flag_deltat_up and not flag_twaiting_up) {
        flag_mc = false;
      } else {
        if (resample_counter >= resample_max) {
          std::cout << "maximum number of resamplings (" << resample_counter
                    << ") reached. stopping..." << std::endl;
          flag_mc = false;
        } else {
          std::cout << "resampling..." << std::endl;
          resample_counter++;

          std::cout << "writing temporary files" << std::endl;
          writeAASequences("temp_" + output_file);
          writeNumericalSequences("temp_" + output_file);
        }
      }
    } else {
      flag_mc = false;
    }
  }

  int idx = output_file.find_last_of(".");
  std::string output_name = output_file.substr(0, idx);
  if (deleteFile("temp_" + output_file) != 0)
    std::cerr << "temporary file deletion failed!" << std::endl;
  else if (deleteFile("temp_" + output_name + "_numerical.txt") != 0)
    std::cerr << "temporary file deletion failed!" << std::endl;
  else if (deleteFile("temp_" + output_name + "_energies.txt") != 0)
    std::cerr << "temporary file deletion failed!" << std::endl;

  std::cout << "writing final sequences... " << std::flush;
  writeAASequences(output_file);
  writeNumericalSequences(output_file);
  std::cout << "done" << std::endl;

  return;
};
