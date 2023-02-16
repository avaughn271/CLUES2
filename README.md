# CLUES

## Purpose
CLUES is a program to infer selection coefficients, evaluate the statistical evidence for selection, and reconstruct historic allele frequencies. CLUES was originally developed by Aaron Stern (see the original paper [here](https://doi.org/10.1371/journal.pgen.1008384)) and is currently maintained by [Andrew Vaughn](https://nielsen-lab.github.io/team/andrew-vaughn) as CLUES2.0. Please report any strange results, errors, or code suggestions to him at [ahv36@berkeley.edu](mailto:ahv36@berkeley.edu).

## Installation

CLUES can be downloaded from this GitHub repo by running the following command:
```bash
$ git clone https://github.com/avaughn271/CLUES2.0
```
CLUES is written in Python and requires the following Python packages:

"numba"

See the following links for installation help:

https://numba.readthedocs.io/en/stable/user/installing.html

### Example commands

TODO

## Input file

CLUES has been tested on two different input files: posterior samples of ARGs (as described in the original CLUES paper) and historic allele samples (as described in TODO). We describe the input files here:

times (TODO)

Samples:

```bash
1.648171304943918258e+01 0.000000000000000000e+00 -inf -inf
160.678571428571 -2.25379492882461 -0.13010868534702 -4.07454193492592
147.035714285714 -0.787457860031187 -0.632993257740198 -4.3428059215206
2.135572624386489338e+01 -inf -inf 0.000000000000000000e+00
...
4.967392959190617603e+02 0.000000000000000000e+00 -inf -inf
4.983198490270729053e+02 0.000000000000000000e+00 -inf -inf
500.0 0.000000000000000000e+00 -inf -inf
```

where the first line is the populations and the subsequent lines are the bi-allelic counts in each population for a number of SNPs. The first and second allele type has no meaning and can be chosen arbitrarily. The population names should only include letters and numbers (no spaces, dashes, underscores, etc.). See the R script "ConvertFromVCF.R" in the "example" folder for a template for converting from VCF files to this input. At minimum, you will need to change the name of the input VCF file and the individual-to-population mapping in this script. Keep in mind that VCF files can be quite complex, and therefore this script may not work for all possible input VCF files. The user should always perform a sanity check between the input and output of this step and should not take the output at face value.

## Running CLUES

CLUES has 3 steps:

`example/example.timeb` is a binary file containing derived/ancestral coalescence times for the SNP of interest. These are estimated using [Relate v1.1](https://myersgroup.github.io/relate/) [[Speidel, *et al.* Nat. Gen. 2019]](https://www.nature.com/articles/s41588-019-0484-x) and processed using the `extract_coals.py` script. 

If you use this program, please cite:
  
  [Stern, *et al.* Plos Gen. (2019)](https://journals.plos.org/plosgenetics/article/metrics?id=10.1371/journal.pgen.1008384) (doi: 10.1371/journal.pgen.1008384)

To find the previous version of `clues`, which uses ARGweaver output (Rasmussen et al, 2014; Hubisz, et al, 2019; [docs here](http://compgen.cshl.edu/ARGweaver/doc/argweaver-d-manual.html)), please go to https://github.com/35ajstern/clues-v0. We are no longer maintaining `clues-v0`
