# ARGinfer

[![Documentation Status](https://readthedocs.org/projects/arginfer/badge/?version=latest)](https://arginfer.readthedocs.io/en/latest/?badge=latest)

``arginfer`` is a probabilistic method for inferring the
Ancestral Recombination Graph (ARG). The method uses a Markov chain
Monte Carlo (MCMC) and samples from the posterior distribution of ARGs
under the coalescent with recombination prior and infinite sites model of mutations.

Please refer to the [documentation](https://arginfer.readthedocs.io/en/latest/) 
for information on installing and using the software.
<!--
# Quick Run
The following command can be used to run 1000 MCMC iterations with
burn-in 200 and retaining every 20 samples (`thining = 20`).
 Also `n` is the number of sequences each `L = 1e5` in length
 evolving in a population of effective size `Ne = 5000`, with
 mutation rate 1e-8 mutations/generation/site and recombination rate
 1e-8 recombinations/generation/site.

```
python3 main.py \
    -I 1000 --thin 20 -b 200 \
    -n 5 -L 1e5 --Ne 5000 \
    -r 1e-8 -mu 1e-8 \
    -O output/out1 \
    --random-seed 1 \
    --plot
```
 The output will be stored in the given path "output/out1". If "--plot"
 is given, the trace plot will be stored in the output path in a "pdf"
 format.
-->
