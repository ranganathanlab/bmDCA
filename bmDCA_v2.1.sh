#!/bin/bash

######################HYPERPARAMETERS SETTINGS

##### BM SETTINGS
LAMBDA_REG1=0.01 #L2 regularization strength for 1p statistics
LAMBDA_REG2=0.01 #L2 regularization strength for 2p statistics
STEP_MAX=2000 #max number of BM steps
ERROR_MAX=0.00001 #exit error
SAVE_PARAMETERS=3 #numero di iterazioni ogni quale si salvano i parametri
STEP_CHECK=$STEP_MAX #numero di iterazioni ogni quale si fa il check

##### LEARNING RATES SETTINGS
EPSILON_0_h=0.01 #starting learning rate for fields
EPSILON_0_J=0.001 #starting learning rate for couplings
ADAPT_UP=1.500 #positive adaptive step for learning rate
ADAPT_DOWN=0.600 # negative adaptive step for learning rate
MIN_STEPH=0.001 # min learning rate for fields
MAX_STEPH=2.50 # max learning rate for fields
MIN_STEPJ=0.00001 # min learning rate for couplings 
MAX_STEPJ_N=2.50 # max learning rate for couplings (to be divided by N)
ERROR_MIN_UPDATE=-1 #minimal number of standard deviations of z varible for having parameter update (if negative or zero all parameters are updated)

##### SAMPLING TIMES SETTINGS
T_WAIT0=10000 # starting thermalization time for MCMC
DELTA_T0=100 # starting sampling time for MCMC
CHECK_ERGO=1 #if 1, routine to check if MC is well thermalized and decorrelated is active
ADAPT_UP_TIME=1.500 #positive adaptive step for sampling/thermalization times
ADAPT_DOWN_TIME=0.600 # negative adaptive step for sampling/thermalization times

##### IMPORTANCE SAMPLING SETTINGS
STEP_IMPORTANCE_MAX=1 #importance sampling maximum iterations 
COHERENCE_MIN=0.9999 #coherence importance sampling

##### MCMC SETTINGS
M=1000 #number of samples for each MCMC run
COUNT_MAX=10 #number of independent MCMC runs

##### PARAMETERS INITIALIZATION SETTINGS
#MC_init=0 #1 se si dispone di un campionamento di equilibio iniziale
#MC_file_input=../PCDSHUFCHECK_REP_M10000_TWF100000_DT1000_MAX0.1MIN0001_L-2_fromPLM_2/MC_analysis_150_long40/MC_samples_ERG.txt
PAR_init=1 #0 choose for user defined parameter initialization, 1 for independent model initialization
Parameters_input=../OUTPUT_smallcoup/parameters_learnt_50.txt #starting parameter file

#LEARNING RATE
LEARN_init=1 #0 choose for user defined learning rates, 1 for homogeneous initialization
LEARN_input=../OUTPUT_test2/learning_rate.txt #starting learning rate file

# SETTINGs FRESH ROUTINE
T_WAIT=$T_WAIT0 # 
DELTA_T=$DELTA_T0 # 
MNEW=$(($M * $COUNT_MAX))
##### CHECK ROUTINE SETTINGS
T_WAIT_CHECK=$T_WAIT # 
DELTA_T_CHECK=$DELTA_T # 
M_CHECK=$M #
COUNT_CHECK=$COUNT_MAX #

##### INPUTS
MSA_file=../$1
echo $MSA_file
Weights_file=../$2
Output_folder=$3
MC_file=MC_samples.txt

#Test_file=../INPUT/msa.txt
#Parameters_true_file=../INPUT/par_l10m2.dat

this_file=bmDCA_v2.1.sh

#####################################################################################

mkdir $Output_folder
cp $this_file $Output_folder/
cd $Output_folder

M_NAT=$(head -1 $MSA_file  | awk '{print $1}')
N=$(head -1 $MSA_file | awk '{print $2}')
Q=$(head -1 $MSA_file | awk '{print $3}')

MAX_STEPJ=$(ps -ef | awk -v MAX_STEPJ_N=$MAX_STEPJ_N -v N=$N 'BEGIN{print MAX_STEPJ_N/N}')

rm out*
rm parameters_learnt*
rm MC_sample*
rm overlap*
rm stat_*
rm ergo*
rm coherence*

error=$ERROR_MAX
step=0

################

#M_TEST=$(head -1 $Test_file  | awk '{print $1}')
#./../compute_energies.out $Test_file  $Parameters_true_file my_energy_true.dat

##### 1 STATISTICS TARGET MSA
 
echo 'Statistics of MSA...'
./../statMSA.out $MSA_file $Weights_file
MEFF=$(awk '{sum+=$1}END{print sum}' $Weights_file)


echo "Sequence length: $N"
echo "Training MSA size: $M_NAT"
echo "Effective Size: $MEFF"
echo "Alphabet size: $Q"
echo "Size of MCMC samples: $MNEW"

##### 2 PARAMETERS INITIALIZATION

if [ $PAR_init == '1' ]
then
echo "Starting point: independent model"

PSEUDOCOUNT_IND=$(ps -ef | awk -v MEFF=$MEFF 'BEGIN{print 1.0/MEFF}')
echo $PSEUDOCOUNT_IND
./../initialize_ind.out  $N $Q parameters_temp.txt stat_align_1p.txt $PSEUDOCOUNT_IND
else
echo "Starting point: user-defined model"
cp $Parameters_input parameters_temp.txt 
fi

if [ $LEARN_init == '1' ]
then
echo "Initializing learning rates: homogeneous"
./../initialize.out  $N $Q learning_rate.txt $EPSILON_0_J $EPSILON_0_h
else
echo "Initializing learning rates: user-defined rates"
cp $LEARN_input learning_rate.txt 
fi
	
./../initialize.out  $N $Q gradient_old.txt 0 0


echo 'stdev_gradh stdev_gradJ gradtot max_gradh max_gradJ stdev_zscore1 stdev_zscore2 zscore_tot %singlepoint_updated %2points_updated sigma_connected_corr correlation_connected_corr slope_connected_corr correlation_single_frequencies' > error.txt
echo '$step, $step_importance, $T_WAIT, $DELTA_T, $auto_corr, $ergo_corr, $cross_corr, $e_start, $e_end, $ergo_err, $auto_err, $e_err, $error, $errorh, $errorj, $DELTA_T, $T_WAIT' > T_wait.txt

##### 3 STARTING BM LOOP

while [ $step -lt $STEP_MAX -a  $error > $ERROR_MAX ]
do
	echo "New Iteration step: $step"
	
	##### 3a-sampling
	flag_mc=0
	while [ $flag_mc == '0' ]
	do
		
		echo $MNEW $N $Q > $MC_file
		#tail -n +2 $MC_file | shuf -n $N_OLD  >> initial_conf.txt
		./../MCMC_rip_v2.out -n $N -q $Q -m $M -T $T_WAIT -t $DELTA_T -s $step -r $COUNT_MAX < parameters_temp.txt 
		tail -n +2 out_samples_montecarlo*.txt  >> $MC_file
				
		rm out_*.txt

		#./../plot_relax.sh
		#mv relax_DE.jpeg relax_DE_$step.jpeg
		
		##### 3b-check correct thermalization and decorrelation

		if [ $CHECK_ERGO == '1' ]
			then	

			./../compute_energies.out $MC_file parameters_temp.txt my_out_energies.txt
			tail -n +1 my_out_energies.txt | awk -v M=$M 'NR%M==1'   > my_energies_start.txt
			tail -n +1 my_out_energies.txt | awk -v M=$M 'NR%M==(M-1)'   > my_energies_end.txt
		
			a=$(awk 'BEGIN{sum=0;sum2=0}{sum+=$1;sum2+=$1*$1;counter++;}END{av=sum/counter; av2=sum2/counter;print counter,av,sqrt(av2-av*av)}' my_energies_start.txt)
			b=$(awk 'BEGIN{sum=0;sum2=0}{sum+=$1;sum2+=$1*$1;counter++;}END{av=sum/counter; av2=sum2/counter;print counter,av,sqrt(av2-av*av)}' my_energies_end.txt)
			echo $a $b  >> my_energies_cfr.txt
			awk -v  COUNT_MAX=$COUNT_MAX '{print $2, $5, sqrt($3*$3+$6*$6)/sqrt(COUNT_MAX)}' my_energies_cfr.txt > my_energies_cfr_err.txt 

			./../autocorrelation.out $MC_file $COUNT_MAX $DELTA_T overlap.txt
			 
			auto_corr=$(tail -1 ergo.txt  | awk {'print $1'})
			check_corr=$(tail -1 ergo.txt  | awk {'print $2'})
			cross_corr=$(tail -1 ergo.txt  | awk {'print $3'})
			cross_check_err=$(tail -1 ergo.txt  | awk {'print $8'})
			auto_cross_err=$(tail -1 ergo.txt  | awk {'print $7'})

			e_start=$(tail -1 my_energies_cfr_err.txt  | awk {'print $1'})
			e_end=$(tail -1 my_energies_cfr_err.txt  | awk {'print $2'})
			e_err=$(tail -1 my_energies_cfr_err.txt  | awk {'print $3'})
		
			echo "ergodicity test:  $check_corr- $cross_corr < $cross_check_err" #if not sat, possible ergodicity issue: sampling times are too short 
			echo "autocorrelation test: $auto_corr - $cross_corr > $auto_cross_err" #if not sat, sampling times are too long
			echo "Twaiting test: $e_start - $e_end < 2*$e_err" #if not sat, waiting time for thermalization is too short
			echo "Twaiting test2: $e_start - $e_end > - 2*$e_err" #if not sat, waiting time for thermalization is too short
			
			flag_deltat_up=$(echo " $check_corr- $cross_corr < $cross_check_err " | bc -l)
			echo $flag_deltat_up
			flag_deltat_down=$(echo " $auto_corr - $cross_corr > $auto_cross_err" | bc -l)
			echo $flag_deltat_down			
			flag_twaiting_up=$(echo " $e_start - $e_end < 2*$e_err" | bc -l)
			echo $flag_twaiting_up
			flag_twaiting_down=$(echo " $e_start - $e_end > -2*$e_err" | bc -l)
			echo $flag_twaiting_down

			if [ $flag_deltat_up == '0' ]
				then
				DELTA_T=$(printf %.0f $(echo "$DELTA_T*$ADAPT_UP_TIME" | bc | tr '.' ','))
				echo "UPDATE DELTAT= $DELTA_T"
			elif  [ $flag_deltat_down == '0' ]	
				then
				DELTA_T=$(printf %.0f $(echo "$DELTA_T*$ADAPT_DOWN_TIME" | bc | tr '.' ',' ))
				echo "UPDATE DELTAT= $DELTA_T"	
			fi

			if [ $flag_twaiting_up == '0' ]
				then
				T_WAIT=$(printf %.0f $(echo "$T_WAIT*$ADAPT_UP_TIME" | bc | tr '.' ','))
				echo "UPDATE T_WAIT=$T_WAIT"
			fi

			if [ $flag_twaiting_down == '0' ]
				then
				T_WAIT=$(printf %.0f $(echo "$T_WAIT*$ADAPT_DOWN_TIME" | bc | tr '.' ','))
				echo "UPDATE T_WAIT=$T_WAIT"
			fi

			if [ $flag_deltat_up == '1' -a $flag_twaiting_up == '1' ]
				then
				flag_mc='1'
			fi

		else
		flag_mc='1'
		fi
	done		

	##### 3c IMPORTANCE SAMPLING LOOP
	step_importance=0
	flag_coherence=1

	cp parameters_temp.txt parameters_ref.txt
	while [ $step_importance -lt $STEP_IMPORTANCE_MAX -a $flag_coherence == '1' ]

	do
			
		let "step_importance = $step_importance+1"
		echo "Importance step: $step_importance"	
		##### 3d-Statistics of MC configuration
		echo 'Statistics of MC...'
		if [ $step_importance -gt 1 ]
			then ./../statMC_sigma_importance_v6.out $MC_file $COUNT_MAX parameters_temp.txt parameters_ref.txt			
			coherence=$(tail -1 coherence_importance.txt  | awk {'print $1'})
			echo $coherence
			flag_coherence=$(echo " $coherence > $COHERENCE_MIN && 1.0/$coherence > $COHERENCE_MIN" | bc -l)
			else  ./../statMC_sigma.out $MC_file $COUNT_MAX parameters_temp.txt
			fi

		#### 3d-estimate likelihood gradient...
		echo 'compute gradient...'
 
		./../compute_error_reparametrization_v3.out stat_MC_1p.txt stat_MC_2p.txt stat_align_1p.txt stat_align_2p.txt $N $Q $ERROR_MAX $LAMBDA_REG1 $LAMBDA_REG2 stat_MC_1p_sigma.txt stat_MC_2p_sigma.txt parameters_temp.txt $ERROR_MIN_UPDATE $MEFF
			
		######check error

		#./../compute_energies.out $Test_file  parameters_temp.txt energy_temp.dat
		#./../compare_energies.out energy_temp.dat my_energy_true.dat $M_test energy_cfr.txt

		##### 3e-update learning rate...
		if [ $step -gt 0 ]	
			then				
			./../update_learning_rate_v4.out $N $Q learning_rate.txt gradient.txt gradient_old.txt $ADAPT_UP $ADAPT_DOWN $MIN_STEPH $MAX_STEPH $MIN_STEPJ $MAX_STEPJ			
			fi		
		cp gradient.txt gradient_old.txt


		##### 3f-compare errors and save

		error=$(tail -1 error.txt | awk {'print $3'})
		errorh=$(tail -1 error.txt | awk {'print $1'})
		errorj=$(tail -1 error.txt | awk {'print $2'})


		echo $step, $step_importance, $T_WAIT, $DELTA_T, $auto_corr, $ergo_corr, $cross_corr, $e_start, $e_end, $ergo_err, $auto_err, $e_err, $error, $errorh, $errorj, $DELTA_T, $T_WAIT >> T_wait.txt
			
		##### 3g check analysis

		let "resto= $step % $STEP_CHECK"
		if [ $resto == '0' -a $step -gt 0 ]
			then
			./../MC_analysis_check.sh $LAMBDA_REG1 $LAMBDA_REG2 $step $T_WAIT_CHECK $DELTA_T_CHECK $M_CHECK $COUNT_CHECK $MSA_file $Weights_file
			#cp MC_analysis_$step/MC_samples_ERG.txt MC_samples.txt
		fi

		##### 3h save parameters
		let "resto= $step % $SAVE_PARAMETERS"
		if [ $resto == '0' ]
			then
			cp parameters_temp.txt parameters_learnt_$step.txt
			cp $MC_file MC_samples_$step.txt
		fi

		let "step = $step+1"
	
		##### 3j-update parameters
		echo 'update parameters...'

		./../update_reparametrization.out $N $Q learning_rate.txt parameters_temp.txt gradient.txt parameters_temp.txt stat_align_1p.txt
		

	done

done

cd ..
