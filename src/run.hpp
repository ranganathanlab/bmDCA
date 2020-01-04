#ifndef BMDCA_RUN_HPP
#define BMDCA_RUN_HPP

#include "mcmc.hpp"
#include "mcmc_stats.hpp"
#include "msa.hpp"
#include "msa_stats.hpp"
#include "utils.hpp"

class Model
{
public:
  potts_model params;
  potts_model learning_rates;
  potts_model gradient;
  int N;
  int Q;

  Model(MSAStats, double, double);

  void writeParams(std::string);
  void writeLearningRates(std::string);
  void writeGradient(std::string);
};

class Sim
{
public:
  Sim(MSAStats msa_stats);
  ~Sim(void);
  void initializeParameters(void);
  void initializeRun(void);
  void computeAutocorrelation(void);
  void checkCurrentStep(void);
  bool computeErrorReparametrization(void);
  void updateLearningRate(void);
  void updateReparameterization(void);
  void run(void);
  void writeData(std::string);
  MSA runMCMC(MSAStats msaStats);

private:
  void readInitialSample(int, int);

  // BM settings
  double lambda_reg1;  // L2 regularization strength for 1p statistics (fields)
  double lambda_reg2;  // L2 regularization strength for 2p statistics (cpling)
  bool use_sca_weight; // whether or not to use rel. entropy for position-
                       // specific regularization
  int step_max;        // max number of BM steps
  double error_max;    // exit error
  int save_parameters; // multiple of iterations at which to save parameters
  int step_check;      // ?

  // Learning parameters
  double epsilon_0_h;      // starting learning rate for fields
  double epsilon_0_J;      // starting learning rate for couplings
  double adapt_up;         // positive adaptive step for learning rate
  double adapt_down;       // negative sdaptive step for learning rate
  double min_step_h;       // min learning rate for fields
  double max_step_h;       // maximum learning rate for fields
  double min_step_J;       // minimum learning rate for couplings
  double max_step_J_N;     // maximum learning rate for couplings (to be
                           // divided by N)
  double error_min_update; // minimal number of standard deviation s of z
                           // variable for having parameter update (if
                           // negative or zero all parameters are updated)

  // Sampling times
  int t_wait_0;           // staring thermalization time for MCMC
  int delta_t_0;          // starging samplign time for MCMC
  bool check_ergo;        // flag to check if MC is well thermalized and
                          // decorrelated
  double adapt_up_time;   // negative adaptive step for sampling/
                          // thermalization times
  double adapt_down_time; // positive adaptive step for sampling/
                          // thermalization times

  // Importance sampling settings
  int step_importance_max; // importance sampling maximum iterations
  double coherence_min;    // coherence importance sampling

  // MCMC settings
  int M;                        // number of samples for each MCMC run
  int count_max;                // number of independent MCMC runs
  bool init_sample = false;     // flag for loading the first positions when
                                // initializing the mcmc from a file
  std::string init_sample_file; // name of file with mcmc initial sample

  // Check routine settings
  int t_wait_check;  // t_wait
  int delta_t_check; // delta_t
  int M_check;       // M
  int count_check;   // count_max

  arma::field<arma::Mat<int>> samples;
  arma::Col<int> initial_sample;

  MSAStats msa_stats;

  // Model model;
  Model* current_model;
  Model* previous_model;

  // MCMC mcmc;
  MCMC* mcmc;
  MCMCStats* mcmc_stats;
};

#endif
