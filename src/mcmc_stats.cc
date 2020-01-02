#include "mcmc_stats.h"

#include <armadillo>
#include <iomanip>
#include <iostream>

#include "utils.h"

// MCMCSamples::MCMCSamples(int reps, int N, int M) {
//   samples = arma::field<arma::Mat<int>>(reps);
//   for (int i = 0; i < reps; i++) {
//     samples(i) = arma::Mat<int>(M, N);
//   }
// }

MCMCStats::MCMCStats(arma::field<arma::Mat<int>> s, potts_model p)
{
  reps = s.n_rows;
  N = s(0).n_cols;
  M = s(0).n_rows;
  Q = 21;
  samples = s;
  params = p;

  computeEnergies();
  // WriteMCMCEnergiesCompat(energies, "my_out_energies.txt");

  computeEnergiesStats();
  // writeEnergyStatsCompat("my_energies_start.txt",
  //                        "my_energies_end.txt",
  //                        "my_energies_cfr.txt",
  //                        "my_energies_cfr_err.txt");
  computeAutocorrelation();
  // writeAutocorrelationStatsCompat("overlap.txt", "overlap_inf.txt",
  // "ergo.txt");
}

void
MCMCStats::updateData(arma::field<arma::Mat<int>> s, potts_model p)
{
  samples = s;
  params = p;

  computeEnergies();
  computeEnergiesStats();
  computeAutocorrelation();
}

void
MCMCStats::computeEnergies(void)
{
  energies = arma::Col<double>(reps * M, arma::fill::zeros);
  double E;
  for (int rep = 0; rep < reps; rep++) {
    for (int seq = 0; seq < M; seq++) {
      E = 0;
      for (int i = 0; i < N; i++) {
        E -= params.h(samples.at(rep)(seq, i), i);
        for (int j = i + 1; j < N; j++) {
          E -=
            params.J.at(i, j)(samples.at(rep)(seq, i), samples.at(rep)(seq, j));
        }
      }
      energies(rep * M + seq) = E;
    }
  }
}

void
MCMCStats::computeEnergiesStats(void)
{
  arma::Col<double> energies_start = arma::Col<double>(reps, arma::fill::zeros);
  arma::Col<double> energies_end = arma::Col<double>(reps, arma::fill::zeros);

  for (int rep = 0; rep < reps; rep++) {
    energies_start.at(rep) = energies(M * rep);
    energies_end.at(rep) = energies(M * rep + M - 1);
  }

  energies_start_avg = arma::mean(energies_start);
  energies_start_sigma = arma::stddev(energies_start, 1); // not sample std dev
  energies_end_avg = arma::mean(energies_end);
  energies_end_sigma = arma::stddev(energies_end, 1); // not sample std dev
  energies_err =
    sqrt((pow(energies_start_sigma, 2) + pow(energies_end_sigma, 2)) / reps);
}

void
MCMCStats::writeEnergyStatsCompat(std::string output_file_start,
                                  std::string output_file_end,
                                  std::string output_file_cfr,
                                  std::string output_file_cfr_err)
{

  std::ofstream output_stream_start(output_file_start);
  output_stream_start << std::fixed << std::setprecision(6);

  std::ofstream output_stream_end(output_file_end);
  output_stream_end << std::fixed << std::setprecision(6);

  std::ofstream output_stream_cfr(output_file_cfr);
  output_stream_cfr << std::fixed << std::setprecision(6);

  std::ofstream output_stream_cfr_err(output_file_cfr_err);
  output_stream_cfr_err << std::fixed << std::setprecision(6);

  for (int rep = 0; rep < reps; rep++) {
    output_stream_start << energies(M * rep) << std::endl;
    output_stream_end << energies(M * rep + M - 1) << std::endl;
  }

  output_stream_cfr << reps << " " << energies_start_avg << " "
                    << energies_start_sigma << " " << reps << " "
                    << energies_end_avg << " " << energies_end_sigma
                    << std::endl;
  output_stream_cfr_err << energies_start_avg << " " << energies_end_avg << " "
                        << energies_err << std::endl;
}

void
MCMCStats::computeAutocorrelation(void)
{
  // overlaps = arma::Col<double>(
  // double *d = (double*)malloc(sizeof(double)*M);
  arma::Col<double> d = arma::Col<double>(M, arma::fill::zeros);
  // double *d2 = (double*)malloc(sizeof(double)*M);
  arma::Col<double> d2 = arma::Col<double>(M, arma::fill::zeros);
  // int *count = (int*)malloc(sizeof(int)*M);
  arma::Col<int> count = arma::Col<int>(M, arma::fill::zeros);
  int id;
  double dinf, dinf2;

  for (int rep = 0; rep < reps; rep++) {
    for (int seq1 = 0; seq1 < M; seq1++) {
      for (int seq2 = seq1 + 1; seq2 < M; seq2++) {
        id = 0;
        for (int i = 0; i < N; i++) {
          if (samples.at(rep)(seq1, i) == samples.at(rep)(seq2, i)) {
            id++;
          }
        }
        d(seq2 - seq1) += (double)id / N;
        d2(seq2 - seq1) += (double)id * id / (N * N);
        count(seq2 - seq1)++;
      }
    }
  }

  dinf = 0;
  dinf2 = 0;
  for (int seq1 = 0; seq1 < M; seq1++) {
    for (int rep1 = 0; rep1 < reps; rep1++) {
      for (int rep2 = rep1 + 1; rep2 < reps; rep2++) {
        id = 0;
        for (int i = 0; i < N; i++) {
          if (samples(rep1)(seq1, i) == samples(rep2)(seq1, i)) {
            id++;
          }
        }
        dinf += (double)id / N;
        dinf2 += (double)(id * id) / (N * N);
      }
    }
  }

  overlaps = arma::Col<double>(M - 2, arma::fill::zeros);
  overlaps_sigma = arma::Col<double>(M - 2, arma::fill::zeros);
  for (int i = 1; i < M - 1; i++) {
    overlaps(i - 1) = (double)d(i) / (double)count(i);
    overlaps_sigma(i - 1) =
      sqrt(1.0 / count(i)) *
      sqrt(d2(i) / (double)(count(i)) - pow(d(i) / (double)(count(i)), 2));
  }

  overlap_inf = 2.0 * dinf / (double)(reps * (reps - 1) * M);
  overlap_inf_sigma =
    sqrt(2.0 / (reps * (reps - 1) * M)) *
    sqrt(2.0 * dinf2 / (double)(reps * (reps - 1) * M) -
         pow(2.0 * dinf / (double)(reps * (reps - 1) * M), 2));

  int i_auto = 1;
  int i_check = maximum((double)M / 10.0, 1.0);

  overlap_cross = (double)2.0 * dinf / (double)(reps * (reps - 1) * M);
  overlap_auto = d(i_auto) / (double)(count(i_auto));
  overlap_check = d(i_check) / (double)(count(i_check));

  sigma_cross = sqrt(2.0 * dinf2 / (double)(reps * (reps - 1) * M) -
                     pow(2.0 * dinf / (double)(reps * (reps - 1) * M), 2));
  sigma_auto = sqrt(d2(i_auto) / (double)(count(i_auto)) -
                    pow(d(i_auto) / (double)(count(i_auto)), 2));
  sigma_check = sqrt(d2(i_check) / (double)(count(i_check)) -
                     pow(d(i_check) / (double)(count(i_check)), 2));

  err_cross_auto = sqrt(pow(sigma_cross, 2) + pow(sigma_auto, 2)) / sqrt(reps);
  err_cross_check =
    sqrt(pow(sigma_cross, 2) + pow(sigma_check, 2)) / sqrt(reps);
  err_check_auto = sqrt(pow(sigma_check, 2) + pow(sigma_auto, 2)) / sqrt(reps);
}

void
MCMCStats::writeAutocorrelationStatsCompat(std::string overlap_file,
                                           std::string overlap_inf_file,
                                           std::string ergo_file)
{
  std::ofstream output_stream_overlap(overlap_file);
  output_stream_overlap << std::fixed << std::setprecision(6);

  std::ofstream output_stream_overlap_inf(overlap_inf_file);
  output_stream_overlap_inf << std::fixed << std::setprecision(6);

  std::ofstream output_stream_ergo(ergo_file);
  output_stream_ergo << std::fixed << std::setprecision(6);

  for (int i = 0; i < M - 2; i++) {
    output_stream_overlap << i << " " << overlaps(i) << " " << overlaps_sigma(i)
                          << std::endl;
  }

  output_stream_overlap_inf << "0 " << overlap_inf << " " << overlap_inf_sigma
                            << std::endl;

  output_stream_ergo << overlap_auto << " " << overlap_check << " "
                     << overlap_cross << " " << sigma_auto << " " << sigma_check
                     << " " << sigma_cross << " " << err_cross_auto << " "
                     << err_cross_check << " " << err_check_auto << std::endl;
};

std::vector<double>
MCMCStats::getEnergiesStats(void)
{
  std::vector<double> values = { energies_start_avg,
                                 energies_start_sigma,
                                 energies_end_avg,
                                 energies_end_sigma,
                                 energies_err };
  return values;
};

std::vector<double>
MCMCStats::getCorrelationsStats(void)
{
  std::vector<double> values = {
    overlap_inf,    overlap_inf_sigma, overlap_auto,  overlap_cross,
    overlap_check,  sigma_auto,        sigma_cross,   sigma_check,
    err_cross_auto, err_cross_check,   err_check_auto
  };
  return values;
};

// arma::Col<double> MCMCStats::getEnergiesStats(void) {
//   arma::Col<double> values = arma::Col<double>(5);
//   values(0) = energies_start_avg;
//   values(1) = energies_start_sigma;
//   values(2) = energies_end_avg;
//   values(3) = energies_end_sigma;
//   values(4) = energies_err;
//   return values;
// };
//
// arma::Col<double> MCMCStats::getCorrelationsStats(void) {
//   arma::Col<double> values = arma::Col<double>(11);
//   values(0) = overlap_inf;
//   values(1) = overlap_inf_sigma;
//   values(2) = overlap_auto;
//   values(3) = overlap_cross;
//   values(4) = overlap_check;
//   values(5) = sigma_auto;
//   values(6) = sigma_cross;
//   values(7) = sigma_check;
//   values(8) = err_cross_auto;
//   values(9) = err_cross_check;
//   values(10) = err_check_auto;
//   return values;
// };

void
MCMCStats::computeSampleStats(void)
{
  frequency_1p = arma::Mat<double>(Q, N, arma::fill::zeros);
  frequency_1p_sigma = arma::Mat<double>(Q, N, arma::fill::zeros);
  frequency_2p = arma::field<arma::Mat<double>>(N, N);
  frequency_2p_sigma = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      frequency_2p.at(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
      frequency_2p_sigma.at(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }

  // for (int pos = 0; pos < M; pos++) {
  //   double* freq_ptr = frequency_1p.colptr(pos);
  //   double* freq_sigma_ptr = frequency_1p_sigma.colptr(pos);
  // }

  arma::Mat<int> n1 = arma::Mat<int>(Q, reps, arma::fill::zeros);
  arma::field<arma::Mat<int>> n2 = arma::field<arma::Mat<int>>(reps);
  for (int rep = 0; rep < reps; rep++) {
    n2.at(rep) = arma::Mat<int>(Q, Q, arma::fill::zeros);
  }

  arma::Col<int> n1av = arma::Col<int>(Q, arma::fill::zeros);
  arma::Mat<int> n2av = arma::Mat<int>(Q, Q, arma::fill::zeros);

  arma::Col<int> n1squared = arma::Col<int>(Q, arma::fill::zeros);
  arma::Mat<int> n2squared = arma::Mat<int>(Q, Q, arma::fill::zeros);

  // output
  // std::cout<<"flag 1"<<std::endl;
  for (int i = 0; i < N; i++) {
    n1.zeros();
    n1av.zeros();
    n1squared.zeros();
    for (int rep = 0; rep < reps; rep++) {
      for (int m = 0; m < M; m++) {
        n1(samples.at(rep).at(m, i), rep)++;
      }
    }
    for (int aa = 0; aa < Q; aa++) {
      for (int rep = 0; rep < reps; rep++) {
        n1av(aa) += n1(aa, rep);
        n1squared(aa) += pow(n1(aa, rep), 2);
      }
      frequency_1p.at(aa, i) = (double)n1av(aa) / M / reps;
      frequency_1p_sigma.at(aa, i) =
        maximum(sqrt(((double)n1squared(aa) / (M * M * reps) -
                      pow((double)n1av(aa) / (M * reps), 2)) /
                     sqrt(reps)),
                0);
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int rep = 0; rep < reps; rep++) {
        n2.at(rep).zeros();
      }
      n2av.zeros();
      n2squared.zeros();
      for (int rep = 0; rep < reps; rep++)
        for (int m = 0; m < M; m++) {
          n2.at(rep)(samples.at(rep).at(m, i), samples.at(rep)(m, j))++;
        }

      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          for (int rep = 0; rep < reps; rep++) {
            n2av.at(aa1, aa2) += n2.at(rep).at(aa1, aa2);
            n2squared.at(aa1, aa2) += pow(n2.at(rep).at(aa1, aa2), 2);
          }
          frequency_2p.at(i, j).at(aa1, aa2) =
            (double)n2av.at(aa1, aa2) / (M * reps);
          frequency_2p_sigma.at(i, j).at(aa1, aa2) =
            maximum(sqrt(((double)n2squared.at(aa1, aa2) / (M * M * reps) -
                          pow((double)n2av.at(aa1, aa2) / (M * reps), 2)) /
                         sqrt(reps)),
                    0);
        }
      }
    }
  }
}

void
MCMCStats::writeFrequency1pCompat(std::string output_file,
                                  std::string output_file_sigma)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);
  std::ofstream output_stream_sigma(output_file_sigma);
  output_stream_sigma << std::fixed << std::setprecision(6);

  for (int i = 0; i < N; i++) {
    output_stream << i;
    output_stream_sigma << i;
    for (int aa = 0; aa < Q; aa++) {
      output_stream << " " << frequency_1p.at(aa, i);
      output_stream_sigma << " " << frequency_1p_sigma.at(aa, i);
    }
    output_stream << std::endl;
    output_stream_sigma << std::endl;
  }
}

void
MCMCStats::writeFrequency2pCompat(std::string output_file,
                                  std::string output_file_sigma)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);
  std::ofstream output_stream_sigma(output_file_sigma);
  output_stream_sigma << std::fixed << std::setprecision(6);

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      output_stream << i << " " << j;
      output_stream_sigma << i << " " << j;
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << " " << frequency_2p.at(i, j).at(aa1, aa2);
          output_stream_sigma << " "
                              << frequency_2p_sigma.at(i, j).at(aa1, aa2);
        }
      }
      output_stream << std::endl;
      output_stream_sigma << std::endl;
    }
  }
}

void
MCMCStats::writeSamplesCompat(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);

  // int reps = samples.n_cols;
  int reps = samples.n_rows;
  int N = samples(0).n_cols;
  int M = samples(0).n_rows;

  output_stream << reps * M << " " << N << " " << 21 << std::endl;

  for (int rep = 0; rep < reps; rep++) {
    for (int i = 0; i < M; i++) {
      output_stream << samples(rep)(i, 0);
      for (int j = 1; j < N; j++) {
        output_stream << " " << samples(rep)(i, j);
      }
      output_stream << std::endl;
    }
  }
};

void
MCMCStats::writeSampleEnergiesCompat(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);
  int M = energies.n_rows;

  for (int i = 0; i < M; i++) {
    output_stream << energies(i) << std::endl;
  }
};
