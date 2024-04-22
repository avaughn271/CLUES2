import numpy, msprime, pyslim, sys

HaploidNe = float(sys.argv[1])
NUMBEROFLEAVESTOSAMPLE = int(sys.argv[2])
sss = float(sys.argv[3])
numtrials = int(sys.argv[4])
L = int(float(sys.argv[5]))
mu = float(sys.argv[6])

if sss == 0.0: sss = 1e-6 # so that the code works properly.

indexx = 1
for filename in range(numtrials):
    while True:
        # define hard sweep model
        sweep_model = msprime.SweepGenicSelection(
            position= 2.0,  # beginning of chromosome
            start_frequency=1.0 / (2 * HaploidNe),
            end_frequency=0.75,
            s=(sss),
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

        IsDerived = [0] * numleaves
        for leaff in range(numleaves):
            currleaf = leaff
            IsDerived[leaff] = currleaf
            while firsttree.parent(currleaf) > -1:
                currleaf =  firsttree.parent(currleaf)
                IsDerived[leaff] = currleaf
        IsDerived = numpy.equal(IsDerived, [max(set(IsDerived), key = IsDerived.count)]) + 0
        mts = msprime.sim_mutations(tsfull, rate = mu, model = "binary", discrete_genome = False) # cahnged this recently
        A = len(mts.tables.mutations)
        B = 0
        for k in mts.variants():
            B = B + 1
        if (len(IsDerived) - numpy.sum(IsDerived)) == round(numleaves/4):
            break

    f = open("Input/TreeTopology" + str(indexx) + ".txt", "w")

    for i in range(len(mts.tables.nodes)):
        f.write( str(i) + "\t" + str(mts.first().parent(i)) + "\n")
    f.close()
    
    with open("Input/Data" + str(indexx) + ".vcf", "w") as vcf_file:
        mts.write_vcf(vcf_file)
    f = open("Input/IsDerived" + str(indexx) + ".txt", "w")
    for i in IsDerived:  f.write(str(i) + "\n")
    f.close()
    
    indexx = indexx + 1
    