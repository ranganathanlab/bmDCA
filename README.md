# bmDCA

Here is the 'working version' code for bmDCA. The main routine
bmDCA\_v2.1.sh contains different scripts (c/c++/awk) and can be
executed in command line in macOS/Linux based operating systems.

Steps to use the code

1/ Compile the code by executing the shell script:

   ./bmDCA\_compile\_v2.1.sh

2/ Pre-process the MSA (fasta format):

   ./bmDCA\_preprocessing.sh [-rw] input\_alignment.fasta

If option -r is used, the input alignment in fasta format will be
converted to a numerical format used by the learning procedure. If
option -w is used, reweighting coefficients will be computed for each
sequence in the alignment. Note that this step may take a long time,
as it is quadratic in the number of sequences of the alignment.
Processed data are then found in the "Processed" folder.

3/ Run bmDCA to learn the model parameters:

  ./bmDCA\_v2.1.sh Processed/msa\_numerical.txt Processed/weights.txt
OutputFolder

The three inputs of bmDCA.sh are:

- Processed/msa\_numerical.txt: this is the target MSA in a numeric
format (see the EXAMPLE folder, the file is generated in the
preprocessing step). The first line contains three integers specifying
the number of sequences, the sequence length and the alphabet size;
- Processed/weights.txt: this is the file containing a single column
with the statistical weights of the MSA sequences, it is generated in
the pre-processing step;
- OutputFolder: this is the folder where all outputs will be saved.

Inside the script bmDCA\_v2.1.sh there are hyperparameters that can be
set, modifying the learning, such as values of regularization or
number of iterations. Inferred parameters are present in the
OutputFolder, in files parameters\_learnt\_%d.txt No stopping procedure
has been implemented to stop the learning. The default number of
iterations of the Boltzmann machine is 2000.

The mapping from amino acids to integers is defined in the following way. Amino acids are ordered as in the following string "-ACDEFGHIKLMNPQRSTVWY". They are then mapped to the integer corresponding to their position in the string, minus one. The gap symbol is mapped to 0, A is mapped to 1, etc ...
