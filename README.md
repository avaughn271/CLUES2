# CLUES2

## Purpose
CLUES2 is a program to infer selection coefficients, evaluate the statistical evidence for selection, and reconstruct historic allele frequencies. It is currently maintained by [Andrew Vaughn](https://nielsen-lab.github.io/team/andrew-vaughn). Please report any strange results, errors, or code suggestions to him at [ahv36@berkeley.edu](mailto:ahv36@berkeley.edu). If you use this program, please cite: [Stern, *et al.* Plos Gen. (2019)](https://doi.org/10.1371/journal.pgen.1008384) and: [Vaughn and Nielsen (2024)](https://academic.oup.com/mbe/article/41/8/msae156/7724092).

## Installation

CLUES2 can be downloaded from this GitHub repo by running the following command:
```bash
$ git clone https://github.com/avaughn271/CLUES2
```
CLUES2 is written in Python and requires the following Python packages:

"numba, scipy, numpy, matplotlib, Bio, pandas"

See the following links for installation help:

https://numba.readthedocs.io/en/stable/user/installing.html

https://scipy.org/install/

https://numpy.org/install/

https://matplotlib.org/stable/index.html

https://biopython.org/wiki/Download

https://pandas.pydata.org/docs/getting_started/install.html

## Important Notes About Parameterization
Population genetics is full of errors and headaches caused by a misplaced factor of 2. In an attempt to avoid any possible confusion by users, I want to clarify 2 possible points about the parameterizations used in CLUES2.

1. We parameterize selection in diploids by setting the relative fitnesses of the 3 genotypes to be 1, 1+hs, and 1+s (h=0.5 by default). This is in concordance with the simulation softwares SLiM and msprime but is different from  the parameterizations of 1, 1+s, 1+2s or  1, 1+s, (1+s)<sup>2</sup> that are often used. This is also different from the parameterization used in the original CLUES method, so do not expect the results to be the same. When analyzing inferred selection coefficients or comparing the results of CLUES2 with other methods, please take note of the parameterizations used and adjust if necessary.

2. In the description of the software here and in the manuscript, N denotes the haploid population size. This is to say that N denotes the number of haplotypes in a given population, NOT the number of diploid individuals. Equivalently, N always denotes the inverse of the coalescence rate. This is in line with the ARG-inference software Relate but is different, for example, from the parameterization used in the simulation software msprime. Therefore, please be aware of the differences in the definitions of N when using CLUES2 with other inference and/or simulation software.

## Running CLUES2

CLUES2 has 3 steps:

## (1) Obtain Input Files

The user will need to process their input data to generate a sample of ancient genotypes and/or a set of samples of coalescence times. The user will need to process ancient genotype samples on their own, but we do  offer conversion scripts from 2 ARG-inference methods: Relate and SINGER. We here briefly describe how to run these conversion scripts, and we include a full description of the input file formats in the "Input File Details" section at the end of this page.

## Relate
If you wish to use Relate to infer ARGs, we provide a script (RelateToCLUES.py) to convert the output of Relate importance sampling for the marginal tree at a given SNP. Please see the documentation of Relate (https://myersgroup.github.io/relate/) [[Speidel, *et al.* Nat. Gen. 2019]](https://www.nature.com/articles/s41588-019-0484-x), specifically the "Sample branch lengths" functionality. Use option "--format n" and set "--first_bp" and "--last_bp" to both be equal to the basepair position of your SNP of interest. We provide a script "RelateToCLUES.py" to convert the output file of this step to the input for CLUES2.

```bash
$ python PATH/RelateToCLUES.py
```

## This step takes as input:

**--RelateSamples** The ".newick" output file of the above Relate step.

**--DerivedFile** A file containing one line for each leaf. Considering the lines as 0-indexed, the i'th line is 0 if leaf i has the ancestral allele and is 1 if leaf i has the derived allele. This can often be read off directly from the VCF file line for the snp of interest or the ".haps" input file to Relate, but be aware that the reference allele and the ancestral allele may be different. The ".mut" file of Relate will describe, for each SNP, which allele it is using as the ancestral allele and whether this SNP was flipped from the input polarization. Please double check this file to make sure that you are using the right allele as the derived allele.

**--out**  The prefix of the output file.

## This step will produce (in the current working directory)

***out*** **_times.txt** A file resembling the input file for derived and ancestral coalescent times, which can be used as input to CLUES2.

The script will also flip the allelic states of the minimum number of leaf nodes as necessary in order to enforce the infinite sites assumption. If no flips are necessary, the message "Infinite sites assumption satisfied. No allele flips necessary." will be printed. Otherwise, we will print the number of total flips we are making and the exact indices of the leaves we are flipping. If allele flips are necessary, we will also produce a file ***out*** **_derived.txt** which will be equivalent to the **--DerivedFile** input file, but with the alleles changed according to the necessary allele flips. If ancient data is used, we also print the number of ancient haplotypes we find with the derived and ancestral alleles. It is always assumed that the youngest haplotype in the tree has age 0, which is to say that the tree represents a mixture of ancient and modern genomes.


## SINGER

CLUES2 can also process the output of SINGER (https://github.com/popgenmethods/SINGER) [[Deng, Nielsen, and Song bioRxiv 2024]](https://www.biorxiv.org/content/10.1101/2024.03.16.585351v1). We provide a script "SingerToCLUES.py" to convert the output of SINGER to the input for CLUES2.

```bash
$ python PATH/SingerToCLUES.py
```

## This step takes as input:

**--position** The genomic position of the mutation you wish to analyze in the SINGER trees.

**--tree_path** The absolute path of the folder containing the sampled SINGER trees. For example, if the user inputs "singeroutput" to this argument, then the script will look for, and calculate coalescence times for,  all of the ".trees" files that are contained in the folder "singeroutput".

**--output**  The prefix of the output file.

## This step will produce (in the current working directory)

***output*** **_times.txt**  A file resembling the input file for derived and ancestral coalescent times, which can be used as input to CLUES2.


## (2) Run Inference

In this step, we run the main CLUES2 program

```bash
$ python PATH/inference.py
```

## This step takes as input:

**--times** The input file of coalescent times as produced, for example, by the above scripts.  You can supply any combination of a **--times** file, an **--ancientSamps** file, and/or an **--ancientHaps** file.   The times in this file, as with all times in CLUES2, should be measured in generations.

**--ancientSamps** The input file of ancient genotype likelihoods as described below. You can supply any combination of a **--times** file, an **--ancientSamps** file, and/or an **--ancientHaps** file. The times in this file, as with all times in CLUES2, should be measured in generations.

**--ancientHaps** The input file of ancient haplotype likelihoods as described below. You can supply any combination of a **--times** file, an **--ancientSamps** file, and/or an **--ancientHaps** file. The times in this file, as with all times in CLUES2, should be measured in generations.

**--popFreq** The modern derived allele frequency. Must be strictly between 0 and 1.

**--N** The effective population size (Ne). This is the HAPLOID effective population size, denoting the number of haplotypes in a given population, NOT the number of diploid individuals. 100 diploid humans corresponds to an N value of 200. Either this or a coal file must be supplied, but not both.

**--coal** The population size file denoting different population sizes through time. Identical to the file format used by (RELATE)[https://myersgroup.github.io/relate/]. The population sizes considered here are the HAPLOID population size, denoting the number of haplotypes in a given population, NOT the number of diploid individuals. 100 diploid humans corresponds to an N value of 200. Either this or a value of N must be supplied, but not both.

**--tCutoff** The maximum time (in generations) to be considered in the analysis.

**--df** (optional) This is the number of discretization points that is used to bin the allele frequencies. A higher number will result in less error due to rounding of allele frequencies but at increased computational cost. We find that having a finer discretization (higher df) is more important when N is large. This is because large population sizes result in smaller allele frequency fluctuations from generation to generation and only a fine discretization grid will be able to model them accurately. The most rigorous way to set df is to steadily increase df until the results appear to have converged, but we find that practically the default value of 450 is sufficient for nearly all cases.

**--timeBins** (optional) A list of epoch breakpoints, sorted in increasing order. These, along with 0 and tCutoff, give the endpoints of the disjoint time intervals in which independent selection coefficients will be inferred. For example, if --timeBins 200 300 is used, 3 separate selection coefficients will be inferred, for each of the intervals [0,200), [200,300), and [300,tCutoff). If this argument is not supplied, one selection coefficient will be inferred for the interval [0,tCutoff). The youngest breakpoint must be greater than 0 and the oldest breakpoint must be less than tCutoff.

**--out** The prefix of the output files.

**--noAlleleTraj** If this flag is used, the inferred allele trajectory will not be estimated. This saves considerable computational time as we do not need to run the Monte Carlo integration and the backward algorithm. Only the inference file will be produced, not the freqs or post files.

**--integration_points**  The number of samples that are used in the  Monte Carlo integration to generate the inferred allele trajectory.  A higher number will result in less variance in the estimate due to finite sampling but at increased computational cost.  Default value is 10 * (number of selection coefficients that are estimated). Example: 10 if no argument is supplied to --timeBins and 30 if  --timeBins 200 300 is used.

**--h**  Dominance coefficient to be used. The relative fitnesses of ancestral allele homozygotes, heterozygotes, and derived allele homozygotes is 1, 1 + hs, 1 + s. Default value is h=0.5, corresponding to additive selection. h is assumed to have the same value through all epochs.


## This step will produce (in the current working directory)

***out_inference.txt*** A file resembling the following:

```bash
logLR -log10(p-value) Epoch1_start Epoch1_end SelectionMLE1
0.3273 0.38 0 1000 0.00162
```
It is a tab-delimited file listing the log-likelihood ratio of selection to no selection, -log10(p) as computed from the log likelihood ratio using the tail of a chi-squared distribution, a list of epoch endpoints, and the MLE estimate of the selection coefficient. If more timeBins are used, the file would look something like this:

```bash
logLR -log10(p-value) Epoch1_start Epoch1_end SelectionMLE1 Epoch2_start Epoch2_end SelectionMLE2 Epoch3_start Epoch3_end SelectionMLE3
176.4062 75.44 0 200 0.01068 200 600 -0.00050 600 5000 -0.00576
```
where selection MLEs are listed to the immediate right of the epoch start and end points.

***out_freqs.txt*** A file containing **df** lines, with each representing one of the discretized allele frequency bins. Not produced if the **--noAlleleTraj**  flag is used.

***out_post.txt*** A file containing **df** lines, with each one containing **tCutoff** comma separated probabilities. The j'th column of the i'th row contains the probability of the derived allele being at frequency freqs[i] at j generations before the present, where freqs represent the frequencies given by the ***out_freqs.txt*** file.  Not produced if the **--noAlleleTraj**  flag is used. 

## (3) Plot Inferred Derived Allele Trajectory.

```bash
$ python PATH/plot_traj.py
```

## This step takes as input:

**--freqs** The input file of frequencies, such as those generated by the previous step (***out_freqs.txt***).

**--post** The input file of posteriors, such as those generated by the previous step (***out_post.txt***).

**--figure** The desired figure prefix.

**--posterior_intervals** (optional) A list of the posterior intervals to plot. Must all be numbers between 0.0 and 1.0, exclusive.  Default value is [0.5, 0.75, 0.95, 0.999].

**--generation_time** (optional) The number of years per generation. Used to set the x-axis labels in units of years. If not specified, the x-axis will be given in generations.

## This step will produce (in the current working directory)

***figure.png***  A figure representing a heatmap of the posterior density of the derived allele frequency. For each generation and each number $x$ in the **--posterior_intervals** list, we plot the interval of frequency bins of minimum width that contains at least $x%$ of the posterior probability. Conceptually, this can be thought of as plotting a set of concentric confidence intervals at each generation. If **--generation_time**  is supplied, then time is converted from generations to years using this conversion factor. Otherwise, time will be plotted in generations.

## Example Commands

We provide 2 example sets of commands for running CLUES2, along with sample input files, in the "examples" folder. 

1. The file AncientHaplotypes, which reproduces the ancient haplotypes analysis of the SLC45A2 locus in the ANA subpopulation, as found in the manuscript:
```bash
python PATHTOCLUES/inference.py --coal example.coal --popFreq  0.98  --ancientHaps  example_haplotypes.csv --out ANA_SLC_output --tCutoff 536 --df 600 --timeBins 89 179
python PATHTOCLUES/plot_traj.py --freqs  ANA_SLC_output_freqs.txt --post ANA_SLC_output_post.txt --figure ANA_SLC_figure  --generation_time 28.0
```
With the first command we run the inference method on the file of ancient haplotypes "example_haplotypes.csv", using coalescence rates as defined by the "example.coal" file. We set "--timeBins 89 179", which means we will infer 3 distinct selection coefficients. With the second command, we plot the inferred derived allele trajectory, assuming a generation time of 28 years per generation.

2. The file ImportanceSampling, which roughly reproduces the importance sampling of gene trees approach used in the pan-ancestry analysis of the MCM6 locus, as found in the manuscript (we say "roughly", because we reduce the number of importance samples used for illustrative purposes):

```bash
PATHTORELATE/scripts/SampleBranchLengths/SampleBranchLengths.sh  --mu 1.25e-8 \
                 -i chrom2 \
                 -o  MCM6 \
                 --coal example.coal
                 --format n \
                 --num_samples 200  \
                 --first_bp  136608646 \
                 --last_bp  136608646 \
                 --distchr2.dist


python PATHTOCLUES/RelateToCLUES.py  --RelateSamples MCM6.newick  --DerivedFile 136608646.txt --out MCM6

python PATHTOCLUES/inference.py --N 30000  --popFreq 0.535 --times MCM6_times.txt  --out ALL_MCM6  --tCutoff 536 --df 600 --noAlleleTraj
```

The first command is our call to the SampleBranchLengths function of Relate, which is used to obtain samples of gene trees at the locus of interest. The second command is where we convert the output of Relate into the input for CLUES2. We take in the output of Relate, along with a file (called  136608646.txt as  136,608,646 is the position of the MCM6 SNP of interest on chromosome 2) which lists which leaves are derived or ancestral. Our RelateToCLUES.py script finds that we need to "flip" 2 leaves (by this we mean changing the allelic state of these leaves) in order to satisfy the infinite sites assumption. Lastly, we run our inference script on the MCM6_times.txt file produced by the previous step. We set haploid Ne to be 30000 and turn on the flag "--noAlleleTraj". This does not produce the files necessary to plot the reconstructed allele frequency trajectory, but speeds up the inference step considerably.

## Input File Details

CLUES2 has been tested on three different input files: posterior samples of ARGs, ancient genotype samples, and ancient haplotype samples. We describe the input files here:

Times: These are the posterior samples of ARGs, more specifically the posterior samples of pairwise coalescence times at the specified SNP. We do not explicitly need the tree topology due to the exchangeability of lineages within the derived and ancestral classes. If you have $M$ samples from the posterior, this file will have $2M$ lines. If we iterate through the samples from  $m=0$ to $m=M-1$, the ($m+1$)'th line will be the coalescence times of derived lineages of the $m$'th sample. The ($m+2$)'th line will be the coalescence times of ancestral lineages of the $m$'th sample. For example, if we sample the following $M=3$ trees from the posterior:


<img src="https://github.com/avaughn271/CLUES2/blob/main/examples/clues1.png">

The file Times would be:

```bash
200,700,850
1000,1250
170,900,1200
50,1300
125,900,2105
83.6,1700
```
We see that there is one entry for each coalescence event in the tree. A coalescence node is a derived node if all of its children have the derived allele **or** one of the immediate child branches of that node is the branch on which the mutation must have arose (marked by a dashed line in the images above). A coalescence event is an ancestral node otherwise. In practice, the algorithm looks at the oldest derived coalescence node and properly treats it as a mixed lineage coalescence event as described in the manuscript, but we label it as a derived node here for simplicity. Note also that the coalescence times in each line must be sorted. You may choose to use only M=1 sample and give CLUES2 a times file with only 2 lines, but you should note that the importance sampling framework will not be used but the given tree will simply be treated as the observed true coalescence tree.

Ancient Times: CLUES2 can also be run on ARGs on ancient data. The input file format is the same as that given above, except for the addition of two lines at the beginning of the file that list the sampling times of the different lineages. In particular, the first line is a semicolon-separated line listing the nonzero sampling times of derived lineages (sampling times of ancient haplotypes only, modern haplotypes are not considered). The second line is a semicolon-separated line listing the nonzero sampling times of ancestral lineages (sampling times of ancient haplotypes only, modern haplotypes are not considered). The subsequent lines are the derived coalescence times of the first importance sample, the ancestral  coalescence times of the first importance sample, the derived coalescence times of the second importance sample, etc. as described above. For example, assume we have 2 modern ancestral lineages, 1 ancestral lineage sampled at time 60, 1 modern derived lineage, and 2 derived lineages sampled at times 40 and 50 respectively. If we sample the following 3 trees (where the parentheses in each leaf node represent the sampling times):

<img src="https://github.com/avaughn271/CLUES2/blob/main/examples/cluesancient1.png">

The file Times would be

```bash
40;50
60
200,700,850
1000,1250
170,900,1200
100,1300
125,900,2105
83.6,1700
```

You read this file into CLUES2 with the --times argument, as when using only modern data. CLUES2 then determines whether ancient samples are being used and correctly parses the input file into the appropriate tree structure.

Ancient Genotype Samples: These input files look like the following file. The first column is the sampling times of the ancient samples (given in generations). The second column is the log genotype likelihood of 0/0 (homozygous ancestral). The third column is the log genotype likelihood of 0/1 (which is to say 1|0 or 0|1). The fourth column is the log genotype likelihood of 1/1 (homozygous derived). For example, the first row of the following file means that an individual was sampled 16.45 generations ago and we are 100% certain has a 0/0 genotype. The second row of the following file means that an individual was sampled 170.6 generations ago to which we assign a likelihood of 0.105 of being 0/0, a likelihood 0.878 of being 0/1, and a likelihood of 0.017 of being 1/1. The data should be sorted in increasing order by sampling time.

```bash
16.45 0.0e+00 -inf -inf
170.6 -2.25379492882461 -0.13010868534702 -4.07454193492592
190.05 -0.787457860031187 -0.632993257740198 -4.3428059215206
213 -inf -inf 0.0
...
496.7 0.0 -inf -inf
498.31 0.0 -inf -inf
500.0 0.0 -inf -inf
```

Ancient Haplotype Samples: These input files look like the following file. The first column is the sampling times of the ancient samples (given in generations). The second column is the log genotype likelihood of the ancestral allele. The third column is the log genotype likelihood of the derived allele. For example, the first row of the following file means that a haplotype was sampled 16.45 generations ago and we are 100% certain that the allele at the SNP of interest was ancestral. The second row of the following file means that an individual was sampled 190.05 generations ago to which we assign a likelihood of 0.1 of being ancestral and a likelihood of 0.9 of being derived.   The data should be sorted in increasing order by sampling time.

```bash
16.45 0.0e+00 -inf
190.05 -2.30258509299 -0.10536051565
...
500.0 0.0 -inf
```

Note that we also allow the genotype likelihoods to be specified in normal (non-log) probabilities as well. This is to say that using the following 2 files are equivalent to using the previous 2 files. 

```bash
16.45 1 0 0
170.6 0.105 0.878 0.017
190.05 0.455 0.531 0.013
213 0 0 1
...
496.7 1 0 0
498.31 1 0 0
500.0 1 0 0
```

```bash
16.45 1 0
190.05 0.1 0.9
...
500.0 1 0
```


