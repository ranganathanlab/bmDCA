#include "model.hpp"

#include <armadillo>

#include "utils.hpp"

Model::Model(std::string parameters_file,
             std::string gradient_file,
             std::string learning_rate_file)
{
  params = loadPottsModelCompat(parameters_file);
  gradient = loadPottsModelCompat(gradient_file);
  learning_rates = loadPottsModelCompat(learning_rate_file);

  N = params.h.n_cols;
  Q = params.h.n_rows;
}

Model::Model(std::string parameters_file_h,
             std::string parameters_file_J,
             std::string gradient_file_h,
             std::string gradient_file_J,
             std::string learning_rate_file_h,
             std::string learning_rate_file_J)
{
  params = loadPottsModel(parameters_file_h, parameters_file_J);
  gradient = loadPottsModel(gradient_file_h, gradient_file_J);
  learning_rates = loadPottsModel(learning_rate_file_h, learning_rate_file_J);

  N = params.h.n_cols;
  Q = params.h.n_rows;
}

Model::Model(MSAStats msa_stats, double epsilon_h, double epsilon_J)
{
  N = msa_stats.getN();
  Q = msa_stats.getQ();
  double pseudocount = 1. / msa_stats.getEffectiveM();

  // Initialize the parameters J and h
  params.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J.at(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
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

  // Initialize the learning rates (epsilon_0_h and epsilon_0_H)
  learning_rates.h = arma::Mat<double>(Q, N);
  learning_rates.h.fill(epsilon_h);
  learning_rates.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      learning_rates.J.at(i, j) = arma::Mat<double>(Q, Q);
      learning_rates.J.at(i, j).fill(epsilon_J);
    }
  }

  // Initialize the gradient
  gradient.J = arma::field<arma::Mat<double>>(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J.at(i, j) = arma::Mat<double>(Q, Q, arma::fill::zeros);
    }
  }
  gradient.h = arma::Mat<double>(Q, N, arma::fill::zeros);
};

void
Model::writeParams(std::string output_file_h, std::string output_file_J)
{
  params.h.save(output_file_h, arma::arma_binary);
  params.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeParamsHDF(std::string output_file)
{
  arma::Mat<int> N_mat = arma::Mat<int>(1, 1);
  arma::Mat<int> Q_mat = arma::Mat<int>(1, 1);
  N_mat.at(0, 0) = N;
  Q_mat.at(0, 0) = Q;
  N_mat.save(arma::hdf5_name(output_file, "N"));
  Q_mat.save(arma::hdf5_name(output_file, "Q", arma::hdf5_opts::append));

  params.h.save(arma::hdf5_name(output_file, "h", arma::hdf5_opts::append));

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      params.J.at(i, j).save(
        arma::hdf5_name(output_file,
                        "J/" + std::to_string(i) + "_" + std::to_string(j),
                        arma::hdf5_opts::append));
    }
  }
};

void
Model::writeParamsCompat(std::string output_file)
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
};

void
Model::writeLearningRates(std::string output_file_h, std::string output_file_J)
{
  learning_rates.h.save(output_file_h, arma::arma_binary);
  learning_rates.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeLearningRatesHDF(std::string output_file)
{
  arma::Mat<int> N_mat = arma::Mat<int>(1, 1);
  arma::Mat<int> Q_mat = arma::Mat<int>(1, 1);
  N_mat.at(0, 0) = N;
  Q_mat.at(0, 0) = Q;
  N_mat.save(arma::hdf5_name(output_file, "N"));
  Q_mat.save(arma::hdf5_name(output_file, "Q", arma::hdf5_opts::append));

  learning_rates.h.save(arma::hdf5_name(output_file, "h", arma::hdf5_opts::append));

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      learning_rates.J.at(i, j).save(
        arma::hdf5_name(output_file,
                        "J/" + std::to_string(i) + "_" + std::to_string(j),
                        arma::hdf5_opts::append));
    }
  }
};

void
Model::writeLearningRatesCompat(std::string output_file)
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
};

void
Model::writeGradient(std::string output_file_h, std::string output_file_J)
{
  gradient.h.save(output_file_h, arma::arma_binary);
  gradient.J.save(output_file_J, arma::arma_binary);
};

void
Model::writeGradientHDF(std::string output_file)
{
  arma::Mat<int> N_mat = arma::Mat<int>(1, 1);
  arma::Mat<int> Q_mat = arma::Mat<int>(1, 1);
  N_mat.at(0, 0) = N;
  Q_mat.at(0, 0) = Q;
  N_mat.save(arma::hdf5_name(output_file, "N"));
  Q_mat.save(arma::hdf5_name(output_file, "Q", arma::hdf5_opts::append));

  gradient.h.save(arma::hdf5_name(output_file, "h", arma::hdf5_opts::append));

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      gradient.J.at(i, j).save(
        arma::hdf5_name(output_file,
                        "J/" + std::to_string(i) + "_" + std::to_string(j),
                        arma::hdf5_opts::append));
    }
  }
};

void
Model::writeGradientCompat(std::string output_file)
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
};
