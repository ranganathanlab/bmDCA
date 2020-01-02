# Boltzmann-machine Direct Coupling Analysis (bmDCA)

**Dependencies:** armadillo c++

## Usage

C implementation of bmDCA adapted from [the
original](https://github.com/matteofigliuzzi/bmDCA) code. Method is described
in:

>  Figliuzzi, M., Barrat-Charlaix, P. & Weigt, M. How Pairwise Coevolutionary
>  Models Capture the Collective Residue Variability in Proteins? Molecular
>  Biology and Evolution 35, 1018–1027 (2018).


Steps to use the code:

### 1. Compile and install the code

To install the program globally (default: `/usr/local`), run:

```
./autogen.sh
make
sudo make install
```

If instead you want to install the code locally, run:
```
./autogen.sh --prefix=${HOME}/.local
make
make install
```

Replace the `--prefix` value with any local path.

In the event you with to uninstall the code, simply run `make uninstall`.

### 2. Pre-process the multiple sequence alignment (FASTA format)

This step is required to convert the MSA text file into numerical format.
```
bmDCA_preprocess.sh -i input_alignment.fasta -d output_directory -r
```

If option `-r` is used, re-weighting coefficients will be computed for each
sequence in the alignment. Note that this step may take a long time, as it is
quadratic in the number of sequences of the alignment. Processed data are then
found in the 'ouput_directory' folder.

### 3. Run bmDCA to learn the model parameters

```
bmDCA_run.sh Processed/msa_numerical.txt Processed/weights.txt OutputFolder
```

The three inputs of bmDCA.sh are:

- *Processed/msa_numerical.txt*: this is the target MSA in a numeric format (see
  the EXAMPLE folder, the file is generated in the preprocessing step). The
  first line contains three integers specifying the number of sequences, the
  sequence length and the alphabet size;
- *Processed/weights.txt*: this is the file containing a single column with the
  statistical weights of the MSA sequences, it is generated in the
  pre-processing step;
- *OutputFolder*: this is the folder where all outputs will be saved.

Inside the script `bmDCA_run.sh` there are hyperparameters that can be set,
modifying the learning, such as values of regularization or number of
iterations. Inferred parameters are present in the `OutputFolder`, in files
`parameters_learnt_%d.txt` No stopping procedure has been implemented to stop
the learning. The default number of iterations of the Boltzmann machine is
2000.

The mapping from amino acids to integers is defined in the following way. Amino
acids are ordered as in the following string "-ACDEFGHIKLMNPQRSTVWY". They are
then mapped to the integer corresponding to their position in the string, minus
one. The gap symbol is mapped to 0, A is mapped to 1, etc...

The output directory contains learned parameters saved every 3 iterations
(default) in files called `parameters_learnt_[it].txt`. Indices of sites in the
sequence go from 0 to L-1 in the output format. The file `error.txt` contained
in the output directory contains information about the fitting quality at
different iterations of the Boltzmann machine, and can be used to decide when
to stop calculations. 

## Example

An example file with processed output is provided in the examples directory. To
use it, run:

```
bmDCA_preprocessing.sh -i example/PF00014_raw.fasta -d processed -r
bmDCA_run.sh processed/msa_numerical.txt processed/weights.txt bminf_example
```

The above commands will first compute sequence weights and convert amino acids
to a numerical code. The results will be stored in the "processed" directory.
The output file should match those within "example/results".

The inference is then run on the two processed files, and the output is stored
in the 'bminf_example' directory.
