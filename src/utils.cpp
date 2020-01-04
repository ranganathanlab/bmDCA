
#include "utils.hpp"

#include <cassert>
#include <string>

SeqRecord::SeqRecord(std::string h, std::string s)
  : header(h)
  , sequence(s){};

void
SeqRecord::print(void)
{
  std::cout << ">" << header << std::endl;
  std::cout << sequence << std::endl;
};

std::string
SeqRecord::getRecord(void)
{
  std::string record_string = ">" + header + "\n" + sequence + "\n";
  return record_string;
};

std::string
SeqRecord::getHeader(void)
{
  return header;
};

std::string
SeqRecord::getSequence(void)
{
  return sequence;
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
