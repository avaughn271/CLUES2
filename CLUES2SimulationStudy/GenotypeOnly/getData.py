
import msprime
import numpy as np
Ne = 10000
L = 1e6  # Length of simulated region

# define hard sweep model
sweep_model = msprime.SweepGenicSelection(
    position=L / 2,  # middle of chrom
    start_frequency=1.0 / (2 * Ne),
    end_frequency=0.5,
    s=0.05,
    dt=1e-6)

reps = msprime.sim_ancestry(
    8,
    model=[sweep_model],#, msprime.StandardCoalescent()],
    population_size=Ne,
    recombination_rate=0,
    sequence_length=L)


print(reps.draw_text())
print(reps.tables)