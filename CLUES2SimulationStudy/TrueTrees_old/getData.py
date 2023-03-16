
import msprime
Ne = 15000

# define hard sweep model
sweep_model = msprime.SweepGenicSelection(
    position= 2.0,  # middle of chrom
    start_frequency=1.0 / (2 * Ne),
    end_frequency=0.5,
    s=0.001,
    dt=1e-6)

reps = msprime.sim_ancestry(
    8,
    model=[sweep_model], #, msprime.StandardCoalescent()
    population_size=Ne,
    recombination_rate=0,
    sequence_length=10)

print(reps.draw_text())
print(reps.tables)