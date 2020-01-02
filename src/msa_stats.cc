#include "msa_stats.h"

#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#define AA_ALPHABET_SIZE 21

MSAStats::MSAStats(MSA msa)
{
  // Initialize
  L = msa.SequenceLength;
  M = msa.SequenceCount;
  Q = AA_ALPHABET_SIZE;

  frequency_1p = arma::Mat<double>(AA_ALPHABET_SIZE, L, arma::fill::zeros);
  frequency_2p = arma::field<arma::Mat<double>>(L, L);
  // frequency_3p = arma::field<arma::Cube<double>>(L, L, L);
  rel_entropy_grad_1p =
    arma::Mat<double>(AA_ALPHABET_SIZE, L, arma::fill::zeros);
  aa_background_frequencies =
    arma::Col<double>::fixed<AA_ALPHABET_SIZE>(arma::fill::zeros);

  aa_background_frequencies = {
    0.000, 0.073, 0.025, 0.050, 0.061, 0.042, 0.072, 0.023, 0.053, 0.064, 0.089,
    0.023, 0.043, 0.052, 0.040, 0.052, 0.073, 0.056, 0.063, 0.013, 0.033
  };
  pseudocount = 0.03;

  // Compute the frequecies (1p statistics) for amino acids (and gaps) for each
  // position. Use pointers to make things speedier.
  M_effective = sum(msa.SequenceWeights);
  int* align_ptr = NULL;
  double* freq_ptr = NULL;
  double* weight_ptr = msa.SequenceWeights.memptr();
  for (int i = 0; i < L; i++) {
    align_ptr = msa.Alignment.colptr(i);
    freq_ptr = frequency_1p.colptr(i);
    for (int m = 0; m < M; m++) {
      *(freq_ptr + *(align_ptr + m)) += *(weight_ptr + m);
    }
  }
  frequency_1p = frequency_1p / M_effective;

  // Compute the 2p statistics
  for (int i = 0; i < L; i++) {
    for (int j = i + 1; j < L; j++) {
      frequency_2p(i, j) = arma::Mat<double>(
        AA_ALPHABET_SIZE, AA_ALPHABET_SIZE, arma::fill::zeros);

      int* align_ptr1 = msa.Alignment.colptr(i);
      int* align_ptr2 = msa.Alignment.colptr(j);
      for (int m = 0; m < M; m++) {
        frequency_2p.at(i, j)(*(align_ptr1 + m), *(align_ptr2 + m)) +=
          *(weight_ptr + m);
      }
      frequency_2p(i, j) = frequency_2p(i, j) / M_effective;
    }
  }

  // // Compute the 3p statistics
  // for (int i = 0; i < L; i++) {
  //   for (int j = i + 1; j < L; j++) {
  //     for (int k = j + 1; k < L; k++) {
  //       frequency_3p(i, j, k) = arma::Cube<double>(
  //         AA_ALPHABET_SIZE, AA_ALPHABET_SIZE, AA_ALPHABET_SIZE,
  //         arma::fill::zeros);
  //
  //       int* align_ptr1 = msa.Alignment.colptr(i);
  //       int* align_ptr2 = msa.Alignment.colptr(j);
  //       int* align_ptr3 = msa.Alignment.colptr(k);
  //       for (int m = 0; m < M; m++) {
  //         frequency_3p.at(i, j, k)(*(align_ptr1 + m), *(align_ptr2 + m),
  //         *(align_ptr3 + m)) +=
  //           *(weight_ptr + m);
  //       }
  //       frequency_3p(i, j, k) = frequency_3p(i, j, k) / M_effective;
  //     }
  //   }
  // }

  std::cout << M << " sequences" << std::endl;
  std::cout << L << " positions" << std::endl;
  std::cout << M_effective << " effective sequences" << std::endl;

  // Update the background frequencies based by computing overall gap frequency
  // theta.
  double theta = 0;
  for (int i = 0; i < L; i++) {
    theta += frequency_1p.at(0, i);
  }
  theta = theta / L;
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
  for (int i = 0; i < L; i++) {
    for (int aa = 0; aa < AA_ALPHABET_SIZE; aa++) {
      pos_freq = tmp(aa, i);
      background_freq = aa_background_frequencies(aa);
      if (pos_freq < 1. && pos_freq > 0.) {
        rel_entropy_grad_1p(aa, i) = log((pos_freq * (1. - background_freq)) /
                                         ((1. - pos_freq) * background_freq));
      }
    }
  }
};

double
MSAStats::GetQ(void)
{
  return Q;
}

double
MSAStats::GetM(void)
{
  return M;
}

double
MSAStats::GetL(void)
{
  return L;
}

double
MSAStats::GetEffectiveM(void)
{
  return M_effective;
}

arma::Mat<double>
MSAStats::GetFrequency1p(void)
{
  return frequency_1p;
}

arma::Mat<double>
MSAStats::GetRelEntropyGradient(void)
{
  return rel_entropy_grad_1p;
}

void
MSAStats::WriteRelEntropyGradientCompat(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);

  for (int i = 0; i < L; i++) {
    for (int aa = 0; aa < AA_ALPHABET_SIZE; aa++) {
      output_stream << i << " " << aa << " " << rel_entropy_grad_1p.at(aa, i)
                    << std::endl;
    }
  }
}

void
MSAStats::WriteFrequency1pCompat(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);

  for (int i = 0; i < L; i++) {
    output_stream << i;
    for (int aa = 0; aa < AA_ALPHABET_SIZE; aa++) {
      output_stream << " " << frequency_1p.at(aa, i);
    }
    output_stream << std::endl;
  }
}

void
MSAStats::WriteFrequency2pCompat(std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);

  for (int i = 0; i < L; i++) {
    for (int j = i + 1; j < L; j++) {
      output_stream << i << " " << j;
      for (int aa1 = 0; aa1 < AA_ALPHABET_SIZE; aa1++) {
        for (int aa2 = 0; aa2 < AA_ALPHABET_SIZE; aa2++) {
          output_stream << " " << frequency_2p(i, j)(aa1, aa2);
        }
      }
      output_stream << std::endl;
    }
  }
}

// void
// MSAStats::WriteFrequency3pCompat(std::string output_file)
// {
//   std::ofstream output_stream(output_file);
//   output_stream << std::fixed << std::setprecision(6);
//
//   for (int i = 0; i < L; i++) {
//     for (int j = i + 1; j < L; j++) {
//       for (int k = j + 1; k < L; k++) {
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
