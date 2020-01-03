#include "msa.hpp"

#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

MSA::MSA(std::string msa_file, bool reweight, double threshold)
{
  readInputMSA(msa_file);
  M = seq_records.size();
  N = getSequenceLength(seq_records.begin()->getSequence());
  makeNumericalMatrix();
  if (reweight) {
    computeSequenceWeights(threshold);
  } else {
    sequence_weights = arma::vec(M, arma::fill::ones);
  }
};

void
MSA::readInputMSA(std::string msa_file)
{
  std::ifstream input_stream(msa_file);

  if (!input_stream) {
    std::cerr << "ERROR: cannot write to '"<< msa_file << "'." << std::endl;
    exit(2);
  }

  /*
   * Read a FASTA-formatted multiple sequence alignment. Each record from the
   * file is stored as a SeqRecord object and appended to the seq_records
   * vector.
   */
  std::string header, sequence, line;
  while (input_stream) {
    input_stream >> line;
    if (line[0] == '>') {
      if (sequence.length() > 0) {
        seq_records.push_back(SeqRecord(header, sequence));
        sequence.clear();
        header.clear();
      }
      header = line;
      header.erase(0, 1);
    } else {
      sequence += line;
      line.clear();
    }
  };
  seq_records.push_back(SeqRecord(header, sequence));
  input_stream.close();
}

void
MSA::makeNumericalMatrix(void)
{
  alignment = arma::Mat<int>(M, N);

  int row_idx = 0;
  for (auto seq = seq_records.begin(); seq != seq_records.end(); seq++) {
    std::string sequence = seq->getSequence();
    int col_idx = 0;
    for (auto aa = sequence.begin(); aa != sequence.end(); aa++) {
      switch (*aa) {
        case '-':
        case 'B':
        case 'J':
        case 'O':
        case 'U':
        case 'X':
        case 'Z':
          alignment(row_idx, col_idx) = 0;
          col_idx++;
          break;
        case 'A':
          alignment(row_idx, col_idx) = 1;
          col_idx++;
          break;
        case 'C':
          alignment(row_idx, col_idx) = 2;
          col_idx++;
          break;
        case 'D':
          alignment(row_idx, col_idx) = 3;
          col_idx++;
          break;
        case 'E':
          alignment(row_idx, col_idx) = 4;
          col_idx++;
          break;
        case 'F':
          alignment(row_idx, col_idx) = 5;
          col_idx++;
          break;
        case 'G':
          alignment(row_idx, col_idx) = 6;
          col_idx++;
          break;
        case 'H':
          alignment(row_idx, col_idx) = 7;
          col_idx++;
          break;
        case 'I':
          alignment(row_idx, col_idx) = 8;
          col_idx++;
          break;
        case 'K':
          alignment(row_idx, col_idx) = 9;
          col_idx++;
          break;
        case 'L':
          alignment(row_idx, col_idx) = 10;
          col_idx++;
          break;
        case 'M':
          alignment(row_idx, col_idx) = 11;
          col_idx++;
          break;
        case 'N':
          alignment(row_idx, col_idx) = 12;
          col_idx++;
          break;
        case 'P':
          alignment(row_idx, col_idx) = 13;
          col_idx++;
          break;
        case 'Q':
          alignment(row_idx, col_idx) = 14;
          col_idx++;
          break;
        case 'R':
          alignment(row_idx, col_idx) = 15;
          col_idx++;
          break;
        case 'S':
          alignment(row_idx, col_idx) = 16;
          col_idx++;
          break;
        case 'T':
          alignment(row_idx, col_idx) = 17;
          col_idx++;
          break;
        case 'V':
          alignment(row_idx, col_idx) = 18;
          col_idx++;
          break;
        case 'W':
          alignment(row_idx, col_idx) = 19;
          col_idx++;
          break;
        case 'Y':
          alignment(row_idx, col_idx) = 20;
          col_idx++;
          break;
      }
    }
    row_idx++;
  }
}

void
MSA::writeMatrix(std::string output_file)
{
  std::ofstream output_stream(output_file);
  // output_stream << seq_records.size() << " " << N << " " << 21
  //               << std::endl;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      if (j + 1 == N) {
        output_stream << alignment(i, j) << std::endl;
      } else {
        output_stream << alignment(i, j) << "\t";
      }
    }
  }
  output_stream << std::endl;
}

void
MSA::printAlignment(void)
{
  for (std::vector<SeqRecord>::iterator it = seq_records.begin();
       it != seq_records.end();
       ++it) {
    it->print();
  }
}

int
MSA::getSequenceLength(std::string sequence)
{
  int valid_aa_count = 0;
  for (std::string::iterator it = sequence.begin(); it != sequence.end();
       ++it) {
    switch (*it) {
      case '-':
      case 'B':
      case 'J':
      case 'O':
      case 'U':
      case 'X':
      case 'Z':
      case 'A':
      case 'C':
      case 'D':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
      case 'K':
      case 'L':
      case 'M':
      case 'N':
      case 'P':
      case 'Q':
      case 'R':
      case 'S':
      case 'T':
      case 'V':
      case 'W':
      case 'Y':
        valid_aa_count += 1;
        break;
    }
  }
  return valid_aa_count;
}

void
MSA::computeSequenceWeights(double threshold)
{
  sequence_weights = arma::vec(M, arma::fill::zeros);
  arma::Mat<int> alignment_T = alignment.t();

  double id;
  int* m1_ptr = nullptr;
  int* m2_ptr = nullptr;
  for (int m1 = 0; m1 < M; ++m1) {
    sequence_weights(m1) += 1;
    m1_ptr = alignment_T.colptr(m1);
    for (int m2 = m1 + 1; m2 < M; ++m2) {
      m2_ptr = alignment_T.colptr(m2);
      id = 0;
      for (int i = 0; i < N; ++i) {
        if (*(m1_ptr + i) == *(m2_ptr + i)) {
          id += 1;
        }
      }
      if (id > threshold * N) {
        sequence_weights(m1) += 1;
        sequence_weights(m2) += 1;
      }
    }
  }

  for (int m1 = 0; m1 < M; ++m1) {
    sequence_weights.at(m1) = 1. / sequence_weights.at(m1);
  }
}

void
MSA::writeSequenceWeights(std::string output_file)
{
  std::ofstream output_stream(output_file);
  for (int i = 0; i < M; i++) {
    output_stream << sequence_weights(i) << std::endl;
  }
}

arma::Mat<int>
GetalignmentMatrix();
