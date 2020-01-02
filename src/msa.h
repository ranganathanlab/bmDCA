#ifndef MSA_H
#define MSA_H

#include <armadillo>
#include <string>
#include <vector>
/* #include <fstream> */

class SeqRecord
{
private:
  const std::string header;
  const std::string sequence;

public:
  SeqRecord(std::string, std::string);
  /* virtual ~SeqRecord(void); */
  void PrintRecord();
  std::string GetRecord();
  std::string GetHeader();
  std::string GetSequence();
};

class MSA
{
public:
  arma::Mat<int> Alignment;
  arma::Col<double> SequenceWeights;
  int SequenceCount;
  int SequenceLength;

  // MSA(std::string);
  MSA(std::string, bool=true, double=0.8);
  void PrintAlignment();
  void WriteMatrixCompat(std::string);
  void WriteSequenceWeightsCompat(std::string);
  /* void SaveAlignment1D(std::string); */
  /* void SaveAlignment2D(std::string); */
  /* virtual ~MSA(void); */

private:
  std::vector<SeqRecord> SeqRecords;
  int GetSequenceLength(std::string);
  void ReadInputMSA(std::string);
  void MakeNumericalMatrix(void);
  void ComputeSequenceWeights(double);
};

#endif
