# Boltzmann-machine Direct Coupling Analysis (bmDCA)

Dependencies (installation instructions detailed below):
 * [Armadillo](http://arma.sourceforge.net/)
 * [GCC] that supports the C++11 standard and
   [OpenMP](https://en.wikipedia.org/wiki/OpenMP)
 * Autotools

This repository contains a C++ reimplementation of bmDCA adapted from [the
original](https://github.com/matteofigliuzzi/bmDCA) code. Method is described
in:

>  Figliuzzi, M., Barrat-Charlaix, P. & Weigt, M. How Pairwise Coevolutionary
>  Models Capture the Collective Residue Variability in Proteins? Molecular
>  Biology and Evolution 35, 1018–1027 (2018).

This code is designed to eliminate the original's excessive file I/O and to
parallelize the MCMC in the inference loop.

## Usage

Steps to use the code:

### 0. Install dependencies

#### Armadillo

Armadillo is a C++ linear algebra library. It can be installed using the
standard package repositories for most Linux distributions (check AUR for Arch
Linux) and Homebrew on macOS (`brew install armadillo`).

__Note for macOS users__: bmDCA depends on `pkg-config` for finding paths for
source files and shared object libraries. The directories where the program
expects to find pkgconfig \*.pc files are listed in the `$PKG_CONFIG_PATH`
variable.  To ensure that pkg-config finds these files for Homebrew-installed
programs, you can append to the variable manually, or you can add the
`pkgconfig_add()` function provided the 'tools/rcparams' file provided in the
repository. Simply append the contents of that file to your shell run
command file (e.g. `.bashrc`), and be sure that it contains the line
`pkgconfig_add /usr/local/Cellar`.

Additionally, the linker can only link libraries found in directories specified
in `LD_LIBRARY_PATH`. To add the armadillo `lib/` directory to this variable,
append to the variable yourself, or you can use the `ld_path_add()` function,
also defined in the 'tools/rcparams' file. Source the function in your shell
run command file and then add `ld_lib_add armadillo`. _This needs to be run
after `pkgconfig_add()`._

__OR__

If a compiled package is not available, you will need to install it
manually. First, make sure that `cmake`, `openblas` (or `blas`), `lapack`,
`arpack`, and `SuperLU` are installed. Then, to download and install armadillo
system wide, run the following (Unix systems only):

```
wget http://sourceforge.net/projects/arma/files/armadillo-9.800.4.tar.xz
tar xf armadillo-9.800.4.tar.xz
cd armadillo-9.800.4
cmake .
make -j4
sudo make install
```

The files will be installed to `/usr/local/include` and `/usr/local/lib` by
default. Make sure that both directories are in your GCC path.

#### GCC

To compile the source code, GCC is recommended. At a minimum, you need a
compiler that supports the C++11 standard. Though optional, one that supports
OpenMP is also recommended. Any GCC later than version 4.2 will suffice.

For Linux users, install a recent GCC from your distributions package
repository.

For Mac users, the default `gcc` is actually to `clang`, which will not allow
compilation with the `-fopenmp` flag. To install a recent version of GCC, run:

```
brew install gcc
```

You have several options to get the `gcc` command to default to the
Homebrew-installed GCC. One option is to alias `gcc` to the path to the new GCC
binary. Code that does this is included in the 'tools/rcparams' file.

#### Autotools

These should already be installed. If not, install `automake` and `pkg-config`
from your system repository.

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

Replace the `--prefix` value with any local path that is part of the system
PATH.

In the event you with to uninstall the code, simply run `sudo make uninstall`
(or `make uninstall` if root permissions aren't needed).

### 2. Run bmDCA to learn the model parameters

This step is required to convert the MSA text file into numerical format.
```
bmdca -i input_alignment.fasta -d output_directory -r -c config_file.conf
```

The command line flags are:
 - `-i`: input MSA, FASTA format
 - `-d`: directory where output is written
 - `-r`: flag to compute re-weighting coefficients for each sequence in the
         alignment (optional)
 - `-c`: config file for bmDCA run hyperparameters (optional)
 - `-t`: threshold for computing default sequence weights (default: `0.8`)
 - `-n`: numerical multiple sequence alignment
 - `-w`: file containing sequence weights

If `-r` is not specified, each sequence will be equally weighted, and if no
config file is supplied, the run will default to hyperparameters hard-coded in
the `initializeParameters()` function defined in `src/run.cpp`. The default
number of iterations for the Boltzmann machine is 2000.

The mapping from amino acids to integers is defined in the following way. Amino
acids are ordered as in the following string "-ACDEFGHIKLMNPQRSTVWY". They are
then mapped to the integer corresponding to their position in the string, minus
one. The gap symbol is mapped to 0, A is mapped to 1, etc...

#### Config file

The fields in the config file:

1. `lambda_reg1` - L2 regularization strength for fields, h (default 0.01)
2. `lambda_reg2` - L2 regularization strength for couplings, J (default 0.01)
3. `step_max` - maximum number of iterations for Boltzmann learning process
   (default 2000)
4. `error_max` - error convergence criterion for stopping (default 1e-05)
5. `save_parameters` - save parameters every `save_parameters` number of steps
   (default 50)
6. `step_check` - _unused_ check learned parameters every `step_check` number
   of steps (default 2000)
7. `epsilon_0_h` - initial learning rate for fields (default 0.01)
8. `epsilon_0_J` - initial learning rate for couplings (default 0.001)
9. `adapt_up` - multiple by which to increase Potts (J and h) gradient (default
   1.5)
10. `adapt_down` - multiple by which to decrease Potts (J and h) gradient
    (default 0.6)
11. `min_step_h` - minimum learning rate for h (default 0.001)
12. `max_step_h` - maximum learning rate for h (default 2.5)
13. `min_step_J` - minimum learning rate for J (default 1e-05)
14. `max_step_J_N` - maximum learning rate for J, scaled by effective number of
    sequences (default 2.5)
15. `error_min_update` - threshold for differences in MSA and MCMC frequencies
    above which parameters (J and h) are updated (default -1)
16. `t_wait_0` - initial burn-in time (default 10000)
17. `delta_t_0` - initial wait time between sampling sequences (default 100)
18. `check_ergo` - flag to check MCMC sample energies and autocorrelations,
    without which wait and burn-in times are not updated (default true)
19. `adapt_up_time` - multiple to increase MCMC wait/burn-in time (default 1.5)
20. `adapt_down_time` - multiple to decrease MCMC wait/burn-in time (default
    0.6)
21. `step_important_max` - maximum number of importance sampling steps
    (default=1 -> importance sampling disabled)
22. `coherence_min` - (default=.9999)
23. `M` - number of sequences to sample for each MCMC replicate (default 1000)
24. `count_max` - number of independent MCMC replicates (default 10)
25. `init_sample` - flag for whether of not to use seed sequence for
    initializing the MCMC (default false)
26. `init_sample_file` - file containing the MCMC seed sequences (default )
27. `temperature` - temperature at which to sample sequences (default 1)
28. `t_wait_check` - _unused_ burn-in time for when running check MCMC (default
    10000)
29. `delta_t_check` - _unused_ wait time between sampled sequences when running
    MCMC check (default 100)
30. `M_check` - number of sequences to sample for each MCMC check replicate
    (default 1000)
31. `count_check` - number of replicates for the MCMC check (default 10)

### 3. Examine the output

The outputs of `bmdca` are:
 - `bmdca_params.conf`: a list of the hyperparameters used in the learning procedure.
 - `energy_%d.dat`: mean and std dev over replicates for sample sequence
   energies at each step of the Markov chain
 - `ergo_%d.dat`: set of autocorrelation calculations for sampled sequences
   used for deciding whether to increase/decrease MCMC wait intervals and
   burn-in times
   1. correlation of sequences 1 wait interval apart.
   2. correlation of sequences M/10 wait intervals apart. (M = # sequences)
   3. cross correlation of sequences
   4. standard deviation of correlations 1 wait intervals apart
   5. standard deviation of correlations M/10 intervals apart
   6. standard deviation of cross correlations
   7. combined deviation of cross and autocorrelations (1 wait interval)
   8. combined deviation of cross and autocorrelations (M/10 wait intervals)
   9. combined deviation of autocorrelations 1 and M/10 intervals apart
 - `MC_energies_%d.txt`: energies of each MCMC sequence, grouped by replicate
 - `MC_samples_%d.txt`: sequences sampled from MCMC, grouped by replicate
 - `msa_numerical.txt`: numerical representation on input MSA
 - `my_energies_cfr_%d.txt`: statistics of energies over replicates, used for
   deciding whether to increase/decrease MCMC wait intervals and burn-in times:
   1. number of replicates
   2. average over replicates of the energies of starting MCMC sequences
   3. standard deviation over replicates of energies of starting MCMC sequences
   4. number of replicates
   5. average over replicates of the energies of ending MCMC sequences
   6. standard deviation over replicates of energies of sending MCMC sequences
 - `my_energies_cfr_err_%d.txt`: additional energies statistics
   1. average over replicates of the energies of starting MCMC sequences
   2. average over replicates of the energies of ending MCMC sequences
   3. combined std dev of energies for staring and ending MCMC sequences
 - `my_energies_end_%d.txt`: energies of ending MCMC sequence for each replicate
 - `my_energies_start_%d.txt`: energies of starting MCMC sequence for each replicate
 - `overlap_%d.txt`: overlap  of pairs of MCMC sequences
   1. number of steps apart (in units of wait time)
   2. mean overlap for all sequences %d steps apart
   3. standard deviation of overlaps for all sequences %d steps apart
 - `parameters_%d.txt`: learned Potts model parameters (J and h)
 - `rel_ent_grad_align_1p.txt`: relative entropy gradient for each amino acid
   at each position
 - `sequence_weights.txt`: weights for each sequence, either a number between 0
   and 1 based on sequence similarity or 1 if re-weighting was not specified
 - `stat_align_1p`: table of frequencies for each amino acid at each position
   in the MSA
 - `stat_align_2p`: table of frequencies for pairs of amino acids at each pair
   of positions in the MSA (due to symmetry, only the 'upper triangle' of
   positions is stored)
 - `stat_MC_1p_%d.txt`: table of frequencies for each amino acid at each
   position of the set of MCMC-sampled sequences.
 - `stat_MC_1p_sigma_%d.txt`: table of standard deviation of frequencies over
   replicates for each amino acid at each position of the set of MCMC-sampled
   sequences.
 - `stat_MC_2p_%d.txt`: table of frequencies for pairs of amino acids at each
   pair of positions from the set of MCMC-sampled sequences
 - `stat_MC_2p_sigma_%d.txt`: table of standard deviation over replicates of
   frequencies for pairs of amino acids at each pair of positions from the set
   of MCMC-sampled sequences

Depending how many times you configure `bmdca` to save steps to disk, the total
data generated can be substantial ( > 1 Gb). At present, the only way to
disable writing of a particular log file is to comment out the code in the
`Sim::run()` function defined in `src/run.cpp`.

Output file formats will probably be changed at a later date.

#### Numerical sequence alignment

This file is a space-delimited file, e.g.:
```
4914 53 21
0 2 10 10 13 16 1 7 6 13 2 1 12 19 17 17 15 19 20 5 18 6 18 18 6 15 2 12 15 5 19 20 6 6 2 7 6 12 9 12 16 5 1 16 4 4 4 2 11 15 18 2 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 5 6 20 9 6 2 6 6 12 9 12 12 5 3 17 3 6 17 2 17 16 10 2 9
```

The first line is:
1. Number of sequences (M)
2. Number of positions (N)
3. Size of amino acid alphabet (all AAs + 1 for gaps) (Q)

#### Learned parameters

The output directory contains learned parameters saved in files called
`parameters_%d.txt`. They contain the parameters for both J and h, formatted
as follows:

```
J [position index i] [position index j] [amino acid index a] [amino acid index b]
.
.
.
h [positoin index i] [amino acid index a]
.
.
.
```

The position indices go from 0 to N-1 (N = # positions), and the amino acid
indices go from 0 to 20 (21 amino acids total, including gaps). 0 corresponds
to a gap.

#### Sequence statistics

The sequence statistics files (e.g. `stat_align_1p.txt` and
`stat_align_2p.txt`) have a different format.

For 1 position (1p) frequencies:
```
[position index] [amino acid frequencies (21)]
.
.
.
```
where `[amino acid frequencies (21)]` is a row of frequencies for each of the
21 positions.

For 2 position (2p) frequencies:
```
[position index i] [position index j] [amino acid frequencies (21x21)]
.
.
.
```
where `[amino acid frequencies (21x21)]` is a row that corresponds to the
frequencies of the `21x21` pairs of amino acids at positions i and j.

## Example

An example FASTA file with processed output is provided in the examples
directory. To use it, run:

```
bmdca -i example/PF00014_raw.fasta -d example/output -r -c example/bmdca.conf
```

A numerical representation and sequence weights of the aforementioned FASTA
file is also available. To use run:

```
bmdca -n example/PF00014_numerical.txt -w example/PF00014_weights.txt \
  -d example/output -c example/bmdca.conf
```
