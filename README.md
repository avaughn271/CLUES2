# CLUES2

## Purpose
CLUES2 is a program to infer selection coefficients, evaluate the statistical evidence for selection, and reconstruct historic allele frequencies. The original CLUES was originally developed by Aaron Stern (see the original paper [here](https://doi.org/10.1371/journal.pgen.1008384)) and is currently maintained by [Andrew Vaughn](https://nielsen-lab.github.io/team/andrew-vaughn) as CLUES2. Please report any strange results, errors, or code suggestions to him at [ahv36@berkeley.edu](mailto:ahv36@berkeley.edu). If you use this program, please cite: [Stern, *et al.* Plos Gen. (2019)](https://doi.org/10.1371/journal.pgen.1008384).

WORK TO DO:
The basic forward backward algorith has been verified. Works well for frequencies, but we now changed the start and end state transitions.

1) Think about allele frequency trajectory, especially considering the boundaries/loss of the mutation
2) Handle the mixed lineage properly
3) Talk about new plotting
4) Evaluate conditional independence of the rest of the arg
5) evaluate importance sampling framework, for example by having relate/argweaver with no recoms and high mutation rate, including checking the proposal density function (why not just comopute same likelihood with s=0. hmmmm.)
6) check how well other arguments work, like the coal file.
7) Then, add the emissions of the args on ancient data, only part that needs to change are the emissions.
8) check coverage probability of posterior intervals.

In terms of improving speed, precompute a lot of the logs and exps. do a better transition probability calculation, as this is actually rather slow. Do Jeff Spence idea 1) bucket means for computing and just shift by mean (or don't shift at all). 2) only compute transition proabilities for square root n number of states, set all else to be 0. Check relative accuracy of this. including for strong selection. Do a better thing for the actual interpolation step. Maybe do a hard index call instead of interpolation. Could even increase number of interpolated points?? Also check what other, slower points of the algorithm are. Also do brent or golden section for search. Check answer is still similar, remember that final function value must also be computed.

Check for a more accurate way of allele transitions, better than diffusion at the boundaries. Maybe diffusion in the middle and something else at the boundaries???

-check coverage probability of the trajectories. Do many independent replicates and pick random point in time. Coverage of all the pooled samples should match up.


## Installation

CLUES can be downloaded from this GitHub repo by running the following command:
```bash
$ git clone https://github.com/avaughn271/CLUES2
```
CLUES is written in Python and requires the following Python packages:

"numba"

See the following links for installation help:

https://numba.readthedocs.io/en/stable/user/installing.html

### Example commands

TODO

## Input file

CLUES has been tested on two different input files: posterior samples of ARGs (as described in the original CLUES paper) and historic allele samples (as described in TODO). We describe the input files here:

Times: These are the posterior samples of ARGs, more specifically the posterior samples of pairwise coalescence times at the specified SNP. We do not explicitly need the tree topology due to the exchangeability of lineages within the derived and ancestral classes. If you have $M$ samples from the posterior, this file will have 2M lines. If we number the samples 0 to $M-1$, the (n+1)th line will be the coalescence times of derived lineages of the nth sample. The (n+2)th line will be the coalescence times of ancestral lineages of the nth sample. For example, if we sample the following M=3 trees from the posterior:

![Test Image 1](example/clues1.png| width=50)
![Test Image 2](example/clues2.png| width=50)
![Test Image 3](example/clues3.png| width=50)
The file Times would look like:

```bash
200,700,850
1000,1250
170,900,1200
50,1300
125,900,2105
83.6,1700
```
We see that there is one entry for each coalescence event in the tree. A coalescence node is a derived node if all of its children have the derived allele **or** one of the immediate child branches of that node is the branch on which the mutation must have arose (marked by a dashed line in the images above). A coalescence event is an ancestral node otherwise. In practice, the algorithm looks at the oldest derived coalescence node and properly treats it as a mixed lineage coalescence event as described in the manuscript, but we label it as a derived node here for simplicity. Note also that the coalescence times in each line must be sorted. You may choose to use only M=1 sample and give CLUES2 a times file with only 2 lines, but you should note that the importance sampling framework will not be used but the given tree will simply be treated as the observed true coalescence tree.


Samples: These input files look like the following file. The first column is the sampling times of the ancient samples (given in generations). The second column is the log genotypes probability of 0/0. The third column is the log genotype probability of 0/1 (which is to say 1|0 or 0|1). The fourth column is the log genotype probability of (1/1). For example, the first row of the following file means that an individual was sampled 16.48 generations ago and we are 100% certain has a 0/0 genotype. The second row of the following file means that an individual was sampled 170.68 generations ago to which we assign a probability of 0.105 of being 0/0, a probability 0.878 of being 0/1, and a probability of 0.017 of being 1/1. Uncertainty in genotype calls of ancient data can be caused by imputation. The data should be sorted in increasing order by sampling time.

```bash
1.648171304943918258e+01 0.000000000000000000e+00 -inf -inf
170.678571428571 -2.25379492882461 -0.13010868534702 -4.07454193492592
190.035714285714 -0.787457860031187 -0.632993257740198 -4.3428059215206
2.135572624386489338e+02 -inf -inf 0.000000000000000000e+00
...
4.967392959190617603e+02 0.000000000000000000e+00 -inf -inf
4.983198490270729053e+02 0.000000000000000000e+00 -inf -inf
500.0 0.000000000000000000e+00 -inf -inf
```

## Running CLUES

CLUES has 3 steps:

## (1) Preprocess Data

The user will need to process their input data to generate a genotype sample or a set of samples of coalescence times. Please see the documentation of Relate [Relate v1.1](https://myersgroup.github.io/relate/) [[Speidel, *et al.* Nat. Gen. 2019]](https://www.nature.com/articles/s41588-019-0484-x) and processed using the `extract_coals.py` script. (TODO)

## (2) Run Inference

In this step, we run the main CLUES program

```bash
$ python PATH/inference.py
```

## This step takes as input:

**--times** The input file of coalescent times as described above. Either this or a sample file must be supplied, but not both.

**--ancientSamps** The input file of ancient genotype proabilities as described above. Either this or a times file must be supplied, but not both.

**--popFreq**  The modern derived allele frequency (corresponding to the allele 1, not 0).

**--N** The population size. This is the HAPLOID population size, denoting the number of haplotypes in a given population, NOT the number of diploid individuals. 100 diploid humans corresponds to an N value of 200. Either this or a coal file must be supplied, but not both.

**--coal** The population size file denoting different population sizes through time. Identical to the file format used by (RELATE)[https://myersgroup.github.io/relate/]. The population sizes considered here are the HAPLOID population size, denoting the number of haplotypes in a given population, NOT the number of diploid individuals. 100 diploid humans corresponds to an N value of 200. Either this or a value of N must be supplied, but not both.

**--tCutoff** The maximum time (in generations) to be considered in the analysis. This is only used if arg samples are used. If ancient allele genotypes are used, the cutoff is automatically set to 1 more than the maximum sampling time.

**--df** This is the number of discretization points that is used to bin the allele frequencies. A higher number will result in less error due to rounding of allele frequencies but at increased computational cost. We find that having a finer discretization (higher df) is more important when N is large and/or s is close to 0. This is because these cases result in smaller allele frequency fluctuations from generation to generation and only a fine discretization grid will be able to model them accurately. The most rigorous way to set df is to steadily increase df until the results appear to have converged, but we find the default value of 400 is sufficient for nearly all cases.

**--timeBins** A file containing one nonnegative number per line, in increasing order. These give the endpoints of the disjoint time intervals in which independent selection coefficients will be inferred. An example file is:

```bash
0.0
50.0
100.0
150.0
```
for which 3 separate selection coefficients will be inferred one each for the intervals [0,50), [50,100), and [100,150). If this file is not supplied, one selection coefficient will be inferred for the interval [0,tCutoff).

**--out** The prefix of the output files.

## This step will produce (in the current working directory)

***out_inference.txt*** A file resembling the following:

```bash
logLR: 122.4629
Epoch	Selection MLE
0-2000	0.00497
```

where the log-likelihood ratio of selection to no selection is listed. The MLE of the selection coefficient is also given for each epoch, as described by the **timeBins** option.

***out_freqs.txt*** A file containing **df** lines, with each representing one of the discretized allele frequency bins. 

***out_post.txt*** A file containing **df** lines, with each representing one containing **tCutoff** comma separated log probailities. The j'th column of the i'th row contains the log-probability of the derived allele being at frequency freqs[i] at j generations before the present, where freqs represent the frequencies given by the ***out_freqs.txt*** file. 

## (3) Plot inferred allele frequency trajectory.

```bash
$ python PATH/plot_traj.py
```

## This step takes as input:

**--freqs** The input file of frequencies, such as those generated by the previous step.

**--post** The input file of log posteriors, such as those generated by the previous step.

**--figure** The desired figure prefix.

## This step will produce (in the current working directory)

***figure.png***  A figure representing a heatmap of the posterior density of derived allele frequency.


