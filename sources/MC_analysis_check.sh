#!/bin/bash

#v11 very high tolerance, flagup -> flag
#v10 change in tw update; change in tolerance
#v9 consistenza tra nomi MSA e parametri
#### MODEL SETTINGS

LAMBDA_REG1=$1 #parametro regolarizzazione 1p
LAMBDA_REG2=$2 #parametro regolarizzazione 2p

#### MONTECARLO SETTINGS
T_WAIT_REL=0 # Twaiting iniziale prima di iniziare il sampling MC
DELTA_T_REL=$5 # distanza tra due configurazioni nel MC
M_REL=$6 #numero di sample per ogni reinizializzazione del MC
COUNT_MAX_REL=$7 #inizializzazioni del MC per ogni step di Boltzmann machine
MC_file_REL=MC_samples_REL.txt
MNEW_REL=$(($M_REL * $COUNT_MAX_REL))

T_WAIT_ERG=$4 # Twaiting iniziale prima di iniziare il sampling MC
DELTA_T_ERG=$5 # distanza tra due configurazioni nel MC
M_ERG=$6 #numero di sample per ogni reinizializzazione del MC
COUNT_MAX_ERG=$7 #inizializzazioni del MC per ogni step di Boltzmann machine
MC_file_ERG=MC_samples_ERG.txt
MNEW_ERG=$(($M_ERG * $COUNT_MAX_ERG))
#Method='MC' #alternatives: MC or PLM 

#PARAMETERS_TRUE_FILE=$8

ERROR_MAX=0.000001
ERROR_MIN_UPDATE=0.000001

SAMPLE_REL=1
SAMPLE_ERG=1

#### INPUTS
folder=MC_analysis_$3
Weights_file=$9
MSA_file=$8
Parameters_input=../parameters_temp.txt

#Test_file=../../PREPROCESSING/msa.txt
Parameters_file=parameters_learnt.txt  




M_NAT=$(head -1 $MSA_file  | awk '{print $1}')
N=$(head -1 $MSA_file | awk '{print $2}')
Q=$(head -1 $MSA_file | awk '{print $3}')
MEFF=$(awk '{sum+=$1}END{print sum}' $Weights_file)
echo $M_NAT


mkdir $folder

cd $folder

####################################

cp $Parameters_input parameters_temp.txt 
cp ../stat_align_1p.txt stat_align_1p.txt
cp ../stat_align_2p.txt stat_align_2p.txt

##############################################

#### Boltzmann  machine

if [ $SAMPLE_ERG == 1 ]
then
	echo $MNEW_ERG $N $Q > $MC_file_ERG
	./../../MCMC_rip_v2.out -n $N -q $Q -m $M_ERG -T $T_WAIT_ERG -t $DELTA_T_ERG -s 0 -r $COUNT_MAX_ERG < parameters_temp.txt
	tail -n +2 out_samples_montecarlo*.txt  >> $MC_file_ERG
	rm out_*.txt
fi

./../../autocorrelation.out $MC_file_ERG $COUNT_MAX_ERG  $DELTA_T_ERG overlap.txt
sh ./../../plot_overlap.sh

###########################

if [ $SAMPLE_REL == 1 ]
then
	echo $MNEW_REL $N $Q > $MC_file_REL
	./../../MCMC_rip_v2.out -n $N -q $Q -m $M_REL -T $T_WAIT_REL -t $DELTA_T_REL -s 0 -r $COUNT_MAX_REL < parameters_temp.txt
	tail -n +2 out_samples_montecarlo*.txt  >> $MC_file_REL
	rm out_*.txt
fi

./../../compute_relaxation.out $MC_file_REL  parameters_temp.txt $COUNT_MAX_REL $DELTA_T_REL energy.dat 		
sh ./../../plot_energy.sh

#########################
		
./../../statMC_sigma.out $MC_file_ERG $COUNT_MAX_ERG parameters_temp.txt			
./../../compute_error_reparametrization_v3.out stat_MC_1p.txt stat_MC_2p.txt stat_align_1p.txt stat_align_2p.txt $N $Q $ERROR_MAX $LAMBDA_REG1 $LAMBDA_REG2 stat_MC_1p_sigma.txt stat_MC_2p_sigma.txt parameters_temp.txt $ERROR_MIN_UPDATE $MEFF
	
echo "stat_MC_1p.txt stat_MC_2p.txt stat_align_1p.txt stat_align_2p.txt $N $Q $ERROR_MAX $LAMBDA_REG1 $LAMBDA_REG2 stat_MC_1p_sigma.txt stat_MC_2p_sigma.txt parameters_temp.txt $ERROR_MIN_UPDATE $MEFF"
sh ./../../plot_stat.sh

#sh ./../../plot_stat_reg.sh

######check error

#./../../compute_energies.out $Test_file  parameters_temp.txt energy_temp_erg.dat
#./../../compute_energies.out $Test_file  $PARAMETERS_TRUE_FILE energy_true_erg.dat
#./../../compare_energies.out energy_temp_erg.dat energy_true_erg.dat $MNEW_ERG energy_cfr_erg.txt


rm parameters_temp.txt
rm stat_*.txt
rm my_corr*

cd ..
