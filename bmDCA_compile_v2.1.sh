#!/bin/bash
echo "Compiling ..."
cd sources
cc reweighting.c -o reweighting.out
cc -o initialize.out initialize.c -lm -O3
cc -o initialize_ind.out initialize_ind.c -lm -O3
cc -o statMSA.out statistics_from_msa.c -lm -O3
cc -o statMC_sigma_importance_v6.out statistics_from_MC_2B_sigma_importance_v6.c -lm -O3 
cc -o statMC_sigma.out statistics_from_MC_2B_sigma.c -lm -O3 
cc -o autocorrelation.out autocorrelation.c -lm -O3 
cc -o compute_error_reparametrization_v3.out compute_error_reparametrization_v3.c -lm -O3 
#cc -o compute_error_v7.out compute_error_v7.c -lm -O3
#cc -o update_v3.out update_v3.c -lm -O3 
cc -o compute_relaxation.out compute_relaxation.c -lm -O3
cc -o update_reparametrization.out update_reparametrization.c -lm -O3
cc -o update_learning_rate_v4.out update_learning_rate_v4.c -lm -O3
#cc -o compare_energies.out compare_energies.c -lm -O3
cc -o compute_energies.out compute_energies.c -lm -O3
g++ MCMC_files_v2/MCMC_rip_v2.cpp  MCMC_files_v2/graph2_rip_init.cpp  MCMC_files_v2/graphs.cpp -o MCMC_rip_v2.out -O3
#chmod +x plot_relax.sh
#chmod +x MC_analysis.sh
chmod +x MC_analysis_check.sh
cd ..
chmod +x bmDCA_preprocessing.sh
chmod +x bmDCA_v2.1.sh
echo "Compilation done."
echo " "
