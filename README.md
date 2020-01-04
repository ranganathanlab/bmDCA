# Boltzmann-machine Direct Coupling Analysis (bmDCA)

Dependencies:
 * [Armadillo](http://arma.sourceforge.net/)
 * [OpenMP](https://en.wikipedia.org/wiki/OpenMP)

## Usage

C++ reimplementation of bmDCA adapted from [the
original](https://github.com/matteofigliuzzi/bmDCA) code. Method is described
in:

>  Figliuzzi, M., Barrat-Charlaix, P. & Weigt, M. How Pairwise Coevolutionary
>  Models Capture the Collective Residue Variability in Proteins? Molecular
>  Biology and Evolution 35, 1018–1027 (2018).

This code is optimized to eliminate the excessive file I/O of the original,
time-consuming printing to std err, and parallelize key steps in the inference
loop.

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

### 2. Run

This step is required to convert the MSA text file into numerical format.
```
bmdca -i input_alignment.fasta -d output_directory -r
```

If option `-r` is used, re-weighting coefficients will be computed for each
sequence in the alignment. Note that this step may take a long time, as it is
quadratic in the number of sequences of the alignment. Processed data are then
found in the 'output\_directory' folder.

### 3. Run bmDCA to learn the model parameters

Inside the source file `src/run.cpp` there are hyperparameters that can be set,
modifying the learning, such as values of regularization or number of
iterations. Inferred parameters are present in the `output_directory`, in files
`parameters_%d.txt` No stopping procedure has been implemented to stop the
learning. The default number of iterations of the Boltzmann machine is 2000.

The mapping from amino acids to integers is defined in the following way. Amino
acids are ordered as in the following string "-ACDEFGHIKLMNPQRSTVWY". They are
then mapped to the integer corresponding to their position in the string, minus
one. The gap symbol is mapped to 0, A is mapped to 1, etc...

The output directory contains learned parameters saved every 3 iterations
(default) in files called `parameters_[it].txt`. Indices of sites in the
sequence go from 0 to L-1 in the output format. The file `error.txt` contained
in the output directory contains information about the fitting quality at
different iterations of the Boltzmann machine, and can be used to decide when
to stop calculations. 

## Example

An example file with processed output is provided in the examples directory. To
use it, run:

```
bmdca -i example/PF00014_raw.fasta -d example/output -r
```
