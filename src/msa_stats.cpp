#include "msa_stats.hpp"

#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>

#define AA_ALPHABET_SIZE 21

MSAStats::MSAStats(MSA msa)
{
  // Initialize
  N = msa.N;
  M = msa.M;
  Q = AA_ALPHABET_SIZE;

  frequency_1p = arma::Mat<double>(AA_ALPHABET_SIZE, N, arma::fill::zeros);
  frequency_2p = arma::field<arma::Mat<double>>(N, N);
  // frequency_3p = arma::field<arma::Cube<double>>(N, N, N);
  rel_entropy_grad_1p =
    arma::Mat<double>(AA_ALPHABET_SIZE, N, arma::fill::zeros);
  aa_background_frequencies =
    arma::Col<double>::fixed<AA_ALPHABET_SIZE>(arma::fill::zeros);

  aa_background_frequencies = {
    0.000, 0.073, 0.025, 0.050, 0.061, 0.042, 0.072, 0.023, 0.053, 0.064, 0.089,
    0.023, 0.043, 0.052, 0.040, 0.052, 0.073, 0.056, 0.063, 0.013, 0.033
  };
  pseudocount = 0.03;

  // Compute the frequecies (1p statistics) for amino acids (and gaps) for each
  // position. Use pointers to make things speedier.
  M_effective = sum(msa.sequence_weights);
  int* align_ptr = nullptr;
  double* freq_ptr = nullptr;
  double* weight_ptr = msa.sequence_weights.memptr();
  for (int i = 0; i < N; i++) {
    align_ptr = msa.alignment.colptr(i);
    freq_ptr = frequency_1p.colptr(i);
    for (int m = 0; m < M; m++) {
      *(freq_ptr + *(align_ptr + m)) += *(weight_ptr + m);
    }
  }
  frequency_1p = frequency_1p / M_effective;

  // Compute the 2p statistics
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      frequency_2p.at(i, j) = arma::Mat<double>(
        AA_ALPHABET_SIZE, AA_ALPHABET_SIZE, arma::fill::zeros);

      int* align_ptr1 = msa.alignment.colptr(i);
      int* align_ptr2 = msa.alignment.colptr(j);
      for (int m = 0; m < M; m++) {
        frequency_2p.at(i, j)(*(align_ptr1 + m), *(align_ptr2 + m)) +=
          *(weight_ptr + m);
      }
      frequency_2p.at(i, j) = frequency_2p.at(i, j) / M_effective;
    }
  }

  // // Compute the 3p statistics
  // for (int i = 0; i < N; i++) {
  //   for (int j = i + 1; j < N; j++) {
  //     for (int k = j + 1; k < N; k++) {
  //       frequency_3p.at(i, j, k) = arma::Cube<double>(
  //         AA_ALPHABET_SIZE, AA_ALPHABET_SIZE, AA_ALPHABET_SIZE,
  //         arma::fill::zeros);
  //
  //       int* align_ptr1 = msa.alignment.colptr(i);
  //       int* align_ptr2 = msa.alignment.colptr(j);
  //       int* align_ptr3 = msa.alignment.colptr(k);
  //       for (int m = 0; m < M; m++) {
  //         frequency_3p.at(i, j, k)(*(align_ptr1 + m), *(align_ptr2 + m),
  //         *(align_ptr3 + m)) +=
  //           *(weight_ptr + m);
  //       }
  //       frequency_3p.at(i, j, k) = frequency_3p.at(i, j, k) / M_effective;
  //     }
  //   }
  // }

  std::cout << M << " sequences" << std::endl;
  std::cout << N << " positions" << std::endl;
  std::cout << Q << " amino acids (including gaps)" << std::endl;
  std::cout << M_effective << " effective sequences" << std::endl;

  // Update the background frequencies based by computing overall gap frequency
  // theta.
  double theta = 0;
  for (int i = 0; i < N; i++) {
    theta += frequency_1p.at(0, i);
  }
  theta = theta / N;
  aa_background_frequencies[0] = theta;
  for (int i = 1; i < AA_ALPHABET_SIZE; i++) {
    aa_background_frequencies[i] = aa_background_frequencies[i] * (1. - theta);
  }

  // Use the positonal and backgrounds frequencies to estimate the relative
  // entropy gradient for each position.
  arma::Mat<double> tmp = frequency_1p * (1. - pseudocount);
  tmp.each_col() += pseudocount * aa_background_frequencies;
  double pos_freq;
  double background_freq;
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < AA_ALPHABET_SIZE; aa++) {
      pos_freq = tmp.at(aa, i);
      background_freq = aa_background_frequencies(aa);
      if (pos_freq < 1. && pos_freq > 0.) {
        rel_entropy_grad_1p.at(aa, i) =
          log((pos_freq * (1. - background_freq)) /
              ((1. - pos_freq) * background_freq));
      }
    }
  }
};

double
MSAStats::getQ(void)
{
  return Q;
};

double
MSAStats::getM(void)
{
  return M;
};

double
MSAStats::getN(void)
{
  return N;
};

double
MSAStats::getEffectiveM(void)
{
  return M_effective;
};

void
MSAStats::writeRelEntropyGradient(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < AA_ALPHABET_SIZE; aa++) {
      output_stream << i << " " << aa << " " << rel_entropy_grad_1p.at(aa, i)
                    << std::endl;
    }
  }
};

void
MSAStats::writeFrequency1p(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int i = 0; i < N; i++) {
    output_stream << i;
    for (int aa = 0; aa < AA_ALPHABET_SIZE; aa++) {
      output_stream << " " << frequency_1p.at(aa, i);
    }
    output_stream << std::endl;
  }
};

void
MSAStats::writeFrequency2p(std::string output_file)
{
  std::ofstream output_stream(output_file);

  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      output_stream << i << " " << j;
      for (int aa1 = 0; aa1 < AA_ALPHABET_SIZE; aa1++) {
        for (int aa2 = 0; aa2 < AA_ALPHABET_SIZE; aa2++) {
          output_stream << " " << frequency_2p(i, j)(aa1, aa2);
        }
      }
      output_stream << std::endl;
    }
  }
};

// void
// MSAStats::writeFrequency3p(std::string output_file)
// {
//   std::ofstream output_stream(output_file);
//
//   for (int i = 0; i < N; i++) {
//     for (int j = i + 1; j < N; j++) {
//       for (int k = j + 1; k < N; k++) {
//         output_stream << i << " " << j << " " << k;
//         for (int aa1 = 0; aa1 < AA_ALPHABET_SIZE; aa1++) {
//           for (int aa2 = 0; aa2 < AA_ALPHABET_SIZE; aa2++) {
//             for (int aa3 = 0; aa3 < AA_ALPHABET_SIZE; aa3++) {
//               output_stream << " " << frequency_3p(i, j, k)(aa1, aa2, aa3);
//             }
//           }
//         }
//         output_stream << std::endl;
//       }
//     }
//   }
// }
