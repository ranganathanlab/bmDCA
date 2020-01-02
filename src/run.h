#ifndef BMDCA_RUN_H
#define BMDCA_RUN_H

#include "mcmc.h"
#include "mcmc_stats.h"
#include "msa.h"
#include "msa_stats.h"
#include "utils.h"

class Model
{
public:
  potts_model params;
  potts_model learning_rates;
  potts_model gradient;

  Model(MSAStats, double, double);
  Model(); // holy fuck is this risky
           // Model(MCMCStats mcmc_stats);
};

class Sim
{
public:
  // Sim(void);
  Sim(MSAStats msa_stats);
  ~Sim(void);
  void initializeParameters(void);
  void initializeRun(void);
  void computeAutocorrelation(void);
  void checkCurrentStep(void);
  // void computeErrorReparametrization(void);
  bool computeErrorReparametrization(void);
  void updateLearningRate(void);
  void updateReparameterization(void);
  void run(void);
  void writeData(int);
  void writeData(void);
  MSA runMCMC(MSAStats msaStats);

private:
  // MCMC *current;
  // MCMC *previous;

  // BM settings
  double lambda_reg1; // L2 regularization strength for 1p statistics
  double lambda_reg2; // L2 regularization strength for 2p statistics
  // double lambda_h;     // field regularization strength
  // double lambda_j;     // couplings regularization strength
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
  int M;         // number of samples for each MCMC run
  int count_max; // number of independent MCMC runs

  //
  // Generated settings
  //
  int t_wait;  // t_wait_0
  int delta_t; // delta_t_0
  int M_new;   // M * count_max

  // Check routine settings
  int t_wait_check;  // t_wait
  int delta_t_check; // delta_t
  int M_check;       // M
  int count_check;   // count_max

  MSAStats msa_stats;

  // Model model;
  Model* current_model;
  Model* previous_model;

  // MCMC mcmc;
  MCMC* mcmc;
  MCMCStats* mcmc_stats;
};

#endif
