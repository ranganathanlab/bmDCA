
#include "utils.h"

#include <iomanip>
#include <iostream>
#include <string>

SeqRecord::SeqRecord(std::string h, std::string s)
  : header(h)
  , sequence(s){};

void
SeqRecord::PrintRecord(void)
{
  std::cout << ">" << header << std::endl;
  std::cout << sequence << std::endl;
};

std::string
SeqRecord::GetRecord(void)
{
  std::string record_string = ">" + header + "\n" + sequence + "\n";
  return record_string;
};

std::string
SeqRecord::GetHeader(void)
{
  return header;
};

std::string
SeqRecord::GetSequence(void)
{
  return sequence;
};

void
WritePottsModelCompat(potts_model model, std::string output_file)
{
  std::ofstream output_stream(output_file);
  output_stream << std::fixed << std::setprecision(6);

  int N = model.h.n_cols;
  int Q = model.h.n_rows;

  // Write J
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      for (int aa1 = 0; aa1 < Q; aa1++) {
        for (int aa2 = 0; aa2 < Q; aa2++) {
          output_stream << "J " << i << " " << j << " " << aa1 << " " << aa2
                        << " " << model.J(i, j)(aa1, aa2) << std::endl;
        }
      }
    }
  }

  // Write h
  for (int i = 0; i < N; i++) {
    for (int aa = 0; aa < Q; aa++) {
      output_stream << "h " << i << " " << aa << " " << model.h(aa, i)
                    << std::endl;
    }
  }
};

double
maximum(double x, double y)
{
  if (x < y)
    return y;
  return x;
};

int
Theta(double x)
{
  if (x > 0)
    return 1;
  return 0;
};

int
Delta(double x)
{
  if (x == 0)
    return 1;
  return 0;
};

double
Max(double a, double b)
{
  if (a > b)
    return a;
  return b;
};

double
Min(double a, double b)
{
  if (a < b)
    return a;
  return b;
};
