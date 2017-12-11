# bmDCA

Here is the 'working version' code for bmDCA. The main routine bmDCA_v2.1.sh is a collage of different scripts (c/c++/awk). 

Steps to use the code
1/ Compiling
./bmDCA_compile_v2.1.sh 

2/ Pre-processing data (if necessary)
./bmDCA_preprocessing.sh [-rw] input_alignment.fasta
If option -r is used, the input alignment in fasta format will be converted to a numerical format used by the learning procedure. 
If option -w is used, reweighting coefficients will be computed for each sequence in the alignment. Note that this step may take a long time, as it is quadratic in the number of sequences of the alignment. 
Processed data are then found in the "Processed" folder. 

3/ Learning the model
./bmDCA_v2.1.sh input_alignment_numerical weights_file output_folder
The three inputs of bmDCA.sh are:
- input_alignment_numerical: this is the target MSA in a numeric format (see the EXAMPLE folder). The first line containts three integers specifying the number of sequences, the sequence length and the alphabet size;
- weights_file: this is the file containing a single column with the statistical weights of the MSA sequences;
- output_folder: is the folder where all outputs will be saved.

Inside the bmDCA_v2.1.sh are hyperparameters that can be set, modifying the learning, such as values of regularizations or number of iterations. 
Inferred parameters are present in the output_folder, in files parameters_learnt_%d.txt
No stopping procedure has been implemented to stop the learning. Default number of iterations of the Boltzmann machine is 2000. 
