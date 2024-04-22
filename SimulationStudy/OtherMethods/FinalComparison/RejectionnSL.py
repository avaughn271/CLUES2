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

def getnewSL(treeseq, isderived):
    numhaps = len(isderived)
    H = numpy.zeros((numhaps, numhaps))
    GenotypeMatrix = treeseq.genotype_matrix()
    Positions = []
    for i in range(len(treeseq.tables.sites)):
        Positions.append(treeseq.tables.sites[i].position)
    for i in range(1,len(Positions)):
        if Positions[i - 1] <= L/2 and Positions[i] >= L/2:
            origleftt = i - 1
            origrightt = i
            break
    for i in range(numhaps - 1):
        for j in range(i + 1,numhaps):
            leftt = origleftt
            rightt = origrightt

            if GenotypeMatrix[rightt,i] != GenotypeMatrix[rightt,j] and GenotypeMatrix[leftt,i] != GenotypeMatrix[leftt,j]:
                H[i,j] = 1
                continue
            elif GenotypeMatrix[rightt,i] != GenotypeMatrix[rightt,j]:
                rightt = rightt - 1
            elif GenotypeMatrix[leftt,i] != GenotypeMatrix[leftt,j]:
                leftt = leftt + 1
            while leftt > 0 and GenotypeMatrix[leftt - 1,i] == GenotypeMatrix[leftt - 1,j]:
                leftt = leftt - 1
            while (rightt < (len(Positions) - 1)) and GenotypeMatrix[rightt + 1,i] == GenotypeMatrix[rightt + 1,j]:
                rightt = rightt + 1
            H[i,j] = rightt - leftt + 2
    nder = numpy.sum(isderived)
    nanc = len(isderived) - numpy.sum(isderived)
    SLA = 0
    SLD = 0
    for i in range(numhaps - 1):
        for j in range(i + 1, numhaps):
            if isderived[i] == 1 and isderived[j] == 1:
                SLD = SLD + H[i,j]
            if isderived[i] == 0 and isderived[j] == 0:
                SLA = SLA + H[i,j]
    SLD = SLD * 2 / (nder * (nder - 1))
    SLA = SLA * 2 / (nanc * (nanc - 1))
    return(numpy.log(SLA/SLD))


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
    IsDerived = [0] * numleaves
    for leaff in range(numleaves):
        currleaf = leaff
        IsDerived[leaff] = currleaf
        while firsttree.parent(currleaf) > -1:
            currleaf =  firsttree.parent(currleaf)
            IsDerived[leaff] = currleaf
    IsDerived = numpy.equal(IsDerived, [max(set(IsDerived), key = IsDerived.count)]) + 0
    if (len(IsDerived) - numpy.sum(IsDerived)) != round(numleaves/4):
        continue
    TajD = getnewSL(mts, IsDerived)

    if abs(TajD  - tajj) < epsilon:
        Totalsels[indexx] = sss
        indexx = indexx + 1
        if indexx == numsamples:
            break
with open("TrueResults/selcoeffs" + str(sys.argv[10]) +".txt", 'w') as fp:
    for item in Totalsels:
        fp.write("%s\n" % item)