#include "run.hpp"

#include <armadillo>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

#include "msa.hpp"
#include "msa_stats.hpp"
#include "utils.hpp"

#define EPSILON 0.00000001

Model::Model(MSAStats msa_stats, double epsilon_h, double epsilon_J)
{
  N = msa_stats.getN();
  Q = msa_stats.getQ();
  double pseudocount = 1. / msa_stats.getEffectiveM();

  params.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J.at(i, j) = arma::Mat<double>(Q, Q);
    }
  }

  params.h = arma::Mat<double>(Q, N, arma::fill::zeros);
  double avg;
  double* freq_ptr = nullptr;
  for (int i = 0; i < N; i++) {
    avg = 0;
    freq_ptr = msa_stats.frequency_1p.colptr(i);
    for (int aa = 0; aa < Q; aa++) {
      avg +=
        log((1. - pseudocount) * (*(freq_ptr + aa)) + pseudocount * (1. / Q));
    }
    for (int aa = 0; aa < Q; aa++) {
      params.h.at(aa, i) =
        log((1. - pseudocount) * (*(freq_ptr + aa)) + pseudocount * (1. / Q)) -
        avg / Q;
    }
  }

  learning_rates.h = arma::Mat<double>(Q, N);
  learning_rates.h.fill(epsilon_h);
  learning_rates.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      learning_rates.J.at(i, j) = arma::Mat<double>(Q, Q);
      learning_rates.J.at(i, j).fill(epsilon_J);
    }
  }

  gradient.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J.at(i, j) = arma::Mat<double>(Q, Q);
    }
  }
  gradient.h = arma::Mat<double>(Q, N, arma::fill::zeros);
};

void
Model::writeParams(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = params.h.n_cols;
  int Q = params.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << params.J.at(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << params.h(aa, i)
                    << std::endl;
    }
  }
}

void
Model::writeLearningRates(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = learning_rates.h.n_cols;
  int Q = learning_rates.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << learning_rates.J.at(i, j)(aa1, aa2)
                        << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << learning_rates.h(aa, i)
                    << std::endl;
    }
  }
}

void
Model::writeGradient(std::string output_file)
{
  std::ofstream output_stream(output_file);

  int N = gradient.h.n_cols;
  int Q = gradient.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << gradient.J.at(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << gradient.h(aa, i)
                    << std::endl;
    }
  }
}

void
Sim::initializeParameters(void)
{
  // BM settings
  lambda_reg1 = 0.01;
  lambda_reg2 = 0.01;
  step_max = 2000;
  error_max = 0.00001;
  save_parameters = 50;
  step_check = step_max;

  // Learning rate settings
  epsilon_0_h = 0.01;
  epsilon_0_J = 0.001;
  adapt_up = 1.5;
  adapt_down = 0.6;
  min_step_h = 0.001;
  max_step_h = 2.5;
  min_step_J = 0.00001;
  max_step_J_N = 2.5; // divide by N later
  error_min_update = -1;

  // sampling time settings
  t_wait_0 = 10000;
  delta_t_0 = 100;
  check_ergo = true;
  adapt_up_time = 1.5;
  adapt_down_time = 0.600;

  // importance sampling settings
  step_importance_max = 1;
  coherence_min = 0.9999;

  // mcmc settings
  M = 1000;            // importance sampling max iterations
  count_max = 10;      // number of independent MCMC runs
  init_sample = false; // flag to load first position for mcmc seqs
  temperature = 1.0;   // temperature at which to sample mcmc

  // check routine settings
  t_wait_check = t_wait_0;
  delta_t_check = delta_t_0;
  M_check = M;
  count_check = count_max;
}

Sim::Sim(MSAStats msa_stats)
  : msa_stats(msa_stats)
{
  initializeParameters();
  current_model = new Model(msa_stats, epsilon_0_h, epsilon_0_J);
  previous_model = new Model(msa_stats, epsilon_0_h, epsilon_0_J);
  mcmc = new MCMC(msa_stats.getN(), msa_stats.getQ());
}

Sim::~Sim(void)
{
  delete current_model;
  delete previous_model;
  delete mcmc;
  delete mcmc_stats;
}

void
Sim::readInitialSample(int N, int Q)
{
  std::ifstream input_stream(init_sample_file);

  if (!input_stream) {
    std::cerr << "ERROR: cannot read '" << init_sample_file << "'."
              << std::endl;
    exit(2);
  }

  std::string line;
  int aa;
  std::getline(input_stream, line);
  for (int n = 0; n < N; n++) {
    std::istringstream iss(line);
    iss >> aa;
    assert(aa < Q);
    initial_sample(n) = aa;
  }
  input_stream.close();
};

void
Sim::run(void)
{

  std::cout << "initializing run" << std::endl << std::endl;

  int N = current_model->N;
  int Q = current_model->Q;

  // Initialize sample data structure
  samples = arma::field<arma::Mat<int>>(count_max);
  for (int i = 0; i < count_max; i++) {
    samples.at(i) = arma::Mat<int>(M, N, arma::fill::zeros);
  }
  mcmc_stats = new MCMCStats(samples, current_model->params);

  if (init_sample) {
    initial_sample = arma::Col<int>(N, arma::fill::zeros);
    readInitialSample(N, Q);
  }

  // BM sampling loop
  int t_wait = t_wait_0;
  int delta_t = delta_t_0;
  for (int step = 0; step <= step_max; step++) {
    std::cout << "Step: " << step << std::endl;

    // Sampling from MCMC (keep trying until correct properties found)
    bool flag_mc = true;
    while (flag_mc) {

      // Draw from MCMC
      std::cout << "loading params to mcmc... ";
      mcmc->load(current_model->params);
      std::cout << "done" << std::endl;

      std::cout << "sampling model with mcmc... " << std::endl;
      if (init_sample) {
        mcmc->sample_init(&samples,
                          count_max,
                          M,
                          N,
                          t_wait,
                          delta_t,
                          &initial_sample,
                          temperature);
      } else {
        mcmc->sample(
          &samples, count_max, M, N, t_wait, delta_t, step, temperature);
      }

      std::cout << "computing mcmc stats... ";
      mcmc_stats->updateData(samples, current_model->params);
      std::cout << "done" << std::endl;

      // Run checks and alter burn-in and wait times
      if (check_ergo) {
        std::vector<double> energy_stats = mcmc_stats->getEnergiesStats();
        std::vector<double> corr_stats = mcmc_stats->getCorrelationsStats();

        double auto_corr = corr_stats.at(2);
        double check_corr = corr_stats.at(3);
        double cross_corr = corr_stats.at(4);
        double cross_check_err = corr_stats.at(9);
        double auto_cross_err = corr_stats.at(8);

        double e_start = energy_stats.at(0);
        double e_end = energy_stats.at(2);
        double e_err = energy_stats.at(4);

        bool flag_deltat_up = false;
        if (check_corr - cross_corr > cross_check_err) {
          flag_deltat_up = true;
        }
        bool flag_deltat_down = false;
        if (auto_corr - cross_corr < auto_cross_err) {
          flag_deltat_down = true;
        }

        bool flag_twaiting_up = false;
        if (e_start - e_end > 2 * e_err) {
          flag_twaiting_up = true;
        }
        bool flag_twaiting_down = false;
        if (e_start - e_end < -2 * e_err) {
          flag_twaiting_down = true;
        }

        if (flag_deltat_up) {
          delta_t = delta_t * adapt_up_time; // truncation: int = double * int;
          std::cout << "increasing wait time to " << delta_t << std::endl;
        } else if (flag_deltat_down) {
          delta_t =
            delta_t * adapt_down_time; // truncation: int = double * int;
          std::cout << "decreasing wait time to " << delta_t << std::endl;
        }

        if (flag_twaiting_up) {
          t_wait = t_wait * adapt_up_time; // truncation: int = double * int;
          std::cout << "increasing burn-in time to " << t_wait << std::endl;
        }
        if (flag_twaiting_down) {
          t_wait = t_wait * adapt_down_time; // truncation: int = double * int;
          std::cout << "decreasing burn-in time to " << t_wait << std::endl;
        }

        if (not flag_deltat_up and not flag_twaiting_up) {
          flag_mc = false;
        } else {
          std::cout << "resampling..." << std::endl;
        }
      } else {
        flag_mc = false;
      }
    }

    // Importance sampling loop
    int step_importance = 0;
    bool flag_coherence = true;
    while (step_importance < step_importance_max and flag_coherence == true) {
      step_importance++;
      if (step_importance > 1) {
        std::cout << "importance sampling step " << step_importance << "... ";
        mcmc_stats->computeSampleStatsImportance(&(current_model->params),
                                                 &(previous_model->params));
        std::cout << "done" << std::endl;

        double coherence = mcmc_stats->Z_ratio;
        if (coherence > coherence_min && 1.0 / coherence > coherence_min) {
          flag_coherence = true;
        } else {
          flag_coherence = false;
        }
      } else {
        std::cout << "computing sequence correlations and energies... ";
        mcmc_stats->computeSampleStats();
        std::cout << "done" << std::endl;
      }

      // Compute error reparametrization
      previous_model->gradient.h = current_model->gradient.h;
      previous_model->gradient.J = current_model->gradient.J;
      std::cout << "computing error and updating gradient... ";
      bool converged = computeErrorReparametrization();
      std::cout << "done" << std::endl;
      if (converged) {
        std::cout << "writing results" << std::endl;
        writeData("final");
        return;
      }

      // Update learning rate
      previous_model->learning_rates.h = current_model->learning_rates.h;
      previous_model->learning_rates.J = current_model->learning_rates.J;
      std::cout << "update learning rate... ";
      updateLearningRate();
      std::cout << "done" << std::endl;

      // Check analysis

      // Save parameters
      if (step % save_parameters == 0 &&
          (step_importance == step_importance_max || flag_coherence == false)) {
        std::cout << "writing step " << step << "... ";
        writeData(std::to_string(step));
        std::cout << "done" << std::endl;
      }

      // Update parameters
      previous_model->params.h = current_model->params.h;
      previous_model->params.J = current_model->params.J;
      std::cout << "updating parameters.. ";
      updateReparameterization();
      std::cout << "done" << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "writing final results... ";
  writeData("final");
  std::cout << "done" << std::endl;
  return;
}

bool
Sim::computeErrorReparametrization(void)
{
  double M_eff = msa_stats.getEffectiveM();
  int N = msa_stats.getN();
  int M = msa_stats.getM();
  int Q = msa_stats.getQ();

  double error_stat_1p = 0;
  double error_stat_2p = 0;
  double error_stat_tot = 0;
  double error_1p = 0;
  double error_2p = 0;
  double error_tot = 0;
  double delta;
  double delta_stat = 0;
  double deltamax_1 = 0;
  double deltamax_2 = 0;
  double rho, beta, den_beta, num_beta, num_rho, den_stat, den_mc, c_mc_av,
    c_stat_av, rho_1p, num_rho_1p, den_stat_1p, den_mc_1p;

  double lambda_h = lambda_reg1;
  double lambda_j = lambda_reg2;

  int count1 = 0;
  int count2 = 0;

  // Compute gradient
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      delta = mcmc_stats->frequency_1p.at(aa, i) -
              msa_stats.frequency_1p.at(aa, i) +
              lambda_h * current_model->params.h.at(aa, i);
      delta_stat =
        (mcmc_stats->frequency_1p.at(aa, i) -
         msa_stats.frequency_1p.at(aa, i)) /
        (pow(msa_stats.frequency_1p.at(aa, i) *
                 (1. - msa_stats.frequency_1p.at(aa, i)) / M_eff +
               pow(mcmc_stats->frequency_1p_sigma.at(aa, i), 2) + EPSILON,
             0.5));
      error_1p += pow(delta, 2);
      error_stat_1p += pow(delta_stat, 2);
      if (pow(delta, 2) > pow(deltamax_1, 2))
        deltamax_1 = sqrt(pow(delta, 2));

      if (sqrt(pow(delta_stat, 2)) > error_min_update) {
        current_model->gradient.h.at(aa, i) = -delta;
        count1++;
      }
    }
  }

  double error_c = 0;
  double c_mc, c_stat = 0;

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          delta = -(msa_stats.frequency_2p.at(i, j).at(aa1, aa2) -
                    mcmc_stats->frequency_2p.at(i, j).at(aa1, aa2) +
                    (mcmc_stats->frequency_1p.at(aa1, i) -
                     msa_stats.frequency_1p.at(aa1, i)) *
                      msa_stats.frequency_1p.at(aa2, j) +
                    (mcmc_stats->frequency_1p.at(aa2, j) -
                     msa_stats.frequency_1p.at(aa2, j)) *
                      msa_stats.frequency_1p.at(aa1, i) -
                    lambda_j * current_model->params.J.at(i, j).at(aa1, aa2));
          delta_stat =
            (mcmc_stats->frequency_2p.at(i, j).at(aa1, aa2) -
             msa_stats.frequency_2p.at(i, j).at(aa1, aa2)) /
            (pow(msa_stats.frequency_2p.at(i, j).at(aa1, aa2) *
                     (1.0 - msa_stats.frequency_2p.at(i, j).at(aa1, aa2)) /
                     M_eff +
                   pow(mcmc_stats->frequency_2p.at(i, j).at(aa1, aa2), 2) +
                   EPSILON,
                 0.5));

          c_mc = mcmc_stats->frequency_2p.at(i, j).at(aa1, aa2) -
                 mcmc_stats->frequency_1p.at(aa1, i) *
                   mcmc_stats->frequency_1p.at(aa2, j);
          c_stat = msa_stats.frequency_2p.at(i, j).at(aa1, aa2) -
                   msa_stats.frequency_1p.at(aa1, i) *
                     msa_stats.frequency_1p.at(aa2, j);
          c_mc_av += c_mc;
          c_stat_av += c_stat;
          error_c += pow(c_mc - c_stat, 2);
          error_2p += pow(delta, 2);
          error_stat_2p += pow(delta_stat, 2);

          if (pow(delta, 2) > pow(deltamax_2, 2)) {
            deltamax_2 = sqrt(pow(delta, 2));
          }
          if (sqrt(pow(delta_stat, 2)) > error_min_update) {
            current_model->gradient.J.at(i, j).at(aa1, aa2) = -delta;
            count2++;
          }
        }
      }
    }
  }

  c_stat_av /= ((N * (N - 1) * Q * Q) / 2);
  c_mc_av /= ((N * (N - 1) * Q * Q) / 2);

  num_rho = num_beta = den_stat = den_mc = den_beta = 0;
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          c_mc = mcmc_stats->frequency_2p.at(i, j).at(aa1, aa2) -
                 mcmc_stats->frequency_1p.at(aa1, i) *
                   mcmc_stats->frequency_1p.at(aa2, j);
          c_stat = msa_stats.frequency_2p.at(i, j).at(aa1, aa2) -
                   msa_stats.frequency_1p.at(aa1, i) *
                     msa_stats.frequency_1p.at(aa2, j);
          num_rho += (c_mc - c_mc_av) * (c_stat - c_stat_av);
          num_beta += (c_mc) * (c_stat);
          den_stat += pow(c_stat - c_stat_av, 2);
          den_mc += pow(c_mc - c_mc_av, 2);
          den_beta += pow(c_stat, 2);
        }
      }
    }
  }

  num_rho_1p = den_stat_1p = den_mc_1p = 0;
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      num_rho_1p += (mcmc_stats->frequency_1p.at(aa, i) - 1.0 / Q) *
                    (msa_stats.frequency_1p.at(aa, i) - 1.0 / Q);
      den_stat_1p += pow(msa_stats.frequency_1p.at(aa, i) - 1.0 / Q, 2);
      den_mc_1p += pow(mcmc_stats->frequency_1p.at(aa, i) - 1.0 / Q, 2);
    }
  }

  beta = num_beta / den_beta;
  rho = num_rho / sqrt(den_mc * den_stat);
  rho_1p = num_rho_1p / sqrt(den_mc_1p * den_stat_1p);

  error_1p = sqrt(error_1p / (N * Q));
  error_2p = sqrt(error_2p / ((N * (N - 1) * Q * Q) / 2));

  error_stat_1p = sqrt(error_stat_1p / (N * Q));
  error_stat_2p = sqrt(error_stat_2p / (N * (N - 1) * Q * Q) / 2);

  error_c = sqrt(error_c / (N * (N - 1) * Q * Q) / 2);

  error_tot = error_1p + error_2p;
  error_stat_tot = error_stat_1p + error_stat_2p;

  bool converged = false;
  if (error_tot < error_max) {
    std::cout << "converged" << std::endl;
    converged = true;
  }
  return converged;
};

void
Sim::updateLearningRate(void)
{
  double M_eff = msa_stats.getEffectiveM();
  int N = msa_stats.getN();
  int M = msa_stats.getM();
  int Q = msa_stats.getQ();
  double max_step_J = max_step_J_N / N;

  double alfa = 0;
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          alfa = Theta(current_model->gradient.J.at(i, j).at(a, b) *
                       previous_model->gradient.J.at(i, j).at(a, b)) *
                   adapt_up +
                 Theta(-current_model->gradient.J.at(i, j).at(a, b) *
                       previous_model->gradient.J.at(i, j).at(a, b)) *
                   adapt_down +
                 Delta(current_model->gradient.J.at(i, j).at(a, b) *
                       previous_model->gradient.J.at(i, j).at(a, b));

          current_model->learning_rates.J.at(i, j).at(a, b) =
            Min(max_step_J,
                Max(min_step_J,
                    alfa * current_model->learning_rates.J.at(i, j).at(a, b)));
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      alfa = Theta(current_model->gradient.h.at(a, i) *
                   previous_model->gradient.h.at(a, i)) *
               adapt_up +
             Theta(-current_model->gradient.h.at(a, i) *
                   previous_model->gradient.h.at(a, i)) *
               adapt_down +
             Delta(current_model->gradient.h.at(a, i) *
                   previous_model->gradient.h.at(a, i));
      current_model->learning_rates.h.at(a, i) =
        Min(max_step_h,
            Max(min_step_h, alfa * current_model->learning_rates.h.at(a, i)));
    }
  }
};

void
Sim::updateReparameterization(void)
{
  double M_eff = msa_stats.getEffectiveM();
  int N = msa_stats.getN();
  int M = msa_stats.getM();
  int Q = msa_stats.getQ();

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int a = 0; a < Q; a++) {
        for (int b = 0; b < Q; b++) {
          current_model->params.J.at(i, j).at(a, b) +=
            current_model->learning_rates.J.at(i, j).at(a, b) *
            current_model->gradient.J.at(i, j).at(a, b);
        }
      }
    }
  };

  arma::Mat<double> Dh = arma::Mat<double>(Q, N, arma::fill::zeros);
  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      for (int j = 0; j < N; j++) {
        if (i < j) {
          for (int b = 0; b < Q; b++) {
            Dh.at(a, i) += -msa_stats.frequency_1p.at(b, j) *
                           current_model->learning_rates.J.at(i, j).at(a, b) *
                           current_model->gradient.J.at(i, j).at(a, b);
          }
        }
        if (i > j) {
          for (int b = 0; b < Q; b++) {
            Dh.at(a, i) += -msa_stats.frequency_1p.at(b, j) *
                           current_model->learning_rates.J.at(j, i).at(b, a) *
                           current_model->gradient.J.at(j, i).at(b, a);
          }
        }
      }
    }
  };

  for (int i = 0; i < N; i++) {
    for (int a = 0; a < Q; a++) {
      current_model->params.h.at(a, i) +=
        current_model->learning_rates.h.at(a, i) *
          current_model->gradient.h.at(a, i) +
        Dh.at(a, i);
    }
  }
};

void
Sim::writeData(std::string id)
{
  current_model->writeParams("parameters_" + id + ".txt");
  mcmc_stats->writeFrequency1p("stat_MC_1p_" + id + ".txt",
                               "stat_MC_1p_sigma_" + id + ".txt");
  mcmc_stats->writeFrequency2p("stat_MC_2p_" + id + ".txt",
                               "stat_MC_2p_sigma_" + id + ".txt");
  mcmc_stats->writeSamples("MC_samples_" + id + ".txt");
  mcmc_stats->writeSampleEnergies("MC_energies_" + id + ".txt");
  mcmc_stats->writeSampleEnergiesRelaxation("energy_" + id + ".dat");
  mcmc_stats->writeEnergyStats("my_energies_start_" + id + ".txt",
                               "my_energies_end_" + id + ".txt",
                               "my_energies_cfr_" + id + ".txt",
                               "my_energies_cfr_err_" + id + ".txt");
  mcmc_stats->writeAutocorrelationStats("overlap_" + id + ".txt",
                                        "overlap_inf_" + id + ".txt",
                                        "ergo_" + id + ".txt");
};
