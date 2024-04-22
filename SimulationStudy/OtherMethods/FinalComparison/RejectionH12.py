import numpy, msprime, pyslim, sys
from numpy import random

HaploidNe = float(sys.argv[1])
NUMBEROFLEAVESTOSAMPLE = int(sys.argv[2])
sss = float(sys.argv[3])
numsamples = int(float(sys.argv[7]))
L = int(float(sys.argv[5]))
mu = float(sys.argv[6])
epsilon = float(sys.argv[9])
tajj = float(sys.argv[8])
print(sys.argv)
if sss == 0.0: sss = 1e-6 # so that the code works properly.

def getH12(treeseq):
    Haplotypes = [[0]]
    bool = True
    for i in range(1,NUMBEROFLEAVESTOSAMPLE):
        for currenthaplotypes in range(len(Haplotypes)):
            if treeseq.diversity([i,Haplotypes[currenthaplotypes][0]]) == 0.0:
                Haplotypes[currenthaplotypes].append(i)
                bool = False
                break
        if bool:
            Haplotypes.append([i])
        bool = True

    Lengths = []
    for i in Haplotypes:
        Lengths.append(len(i))
    Lengths.sort()
    tot = 0
    for i in Lengths:
        tot = tot + i/sum(Lengths) * i/sum(Lengths)
    tot = tot + 2 * Lengths[len(Lengths) - 1]/sum(Lengths) *  Lengths[len(Lengths) - 2]/sum(Lengths)
    return(tot)

Totalsels = [-1] * numsamples
indexx = 0
while True:
    sss = random.uniform(0, 0.025)
    # define hard sweep model
    sweep_model = msprime.SweepGenicSelection(
        position= 2.0,  # beginning of chromosome
        start_frequency=1.0 / (2 * HaploidNe),
        end_frequency=0.75,
        s=(sss*2),
        dt=1e-6)

    reps = msprime.sim_ancestry(
        NUMBEROFLEAVESTOSAMPLE/2,
        model=[sweep_model],
        population_size=HaploidNe/2.0,
        recombination_rate=0.0,
        sequence_length=L)

    demm = msprime.Demography.from_tree_sequence(reps)
    for pop in demm.populations:
        pop.initial_size = HaploidNe / 2.0  ##POSITON1 This is diploid size

    tsfull = pyslim.recapitate(reps, demography= demm, recombination_rate=0.0)

    firsttree = reps.first()

    numleaves = 0
    for k in firsttree.leaves():
        numleaves = numleaves + 1

    mts = msprime.sim_mutations(tsfull, rate = mu, model = "binary", discrete_genome = False) # cahnged this recently
    TajD = getH12(mts)
    if abs(TajD  - tajj) < epsilon:
        Totalsels[indexx] = sss
        indexx = indexx + 1
        if indexx == numsamples:
            break
with open("TrueResults/selcoeffs" + str(sys.argv[10]) +".txt", 'w') as fp:
    for item in Totalsels:
        fp.write("%s\n" % item)