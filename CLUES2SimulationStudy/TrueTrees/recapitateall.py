import numpy, tskit, msprime, pyslim, warnings, os
NUMBEROFLEAVESTOSAMPLE = 100
#If there are 10 nodes, there are 9 total coalescences. However, we do not track the
#first mixed lineage. This gives 8 total coalescence times.
from random import sample
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
files = os.listdir("SLIMTREES")
indexx = 1
for filename in files:
    tstemp = tskit.load("./SLIMTREES/" + filename)

    demm = msprime.Demography.from_tree_sequence(tstemp)
    for pop in demm.populations:
        pop.initial_size = 15000.0 ##POSITON1 This is diploid size

    tsfull = pyslim.recapitate(tstemp, demography= demm, #weird results with ancestral Ne = 15000
                            recombination_rate=0.0) ##POSITION2 RANDOM SEED

    firsttree = (tsfull.first()).split_polytomies()
    orig_max_roots = max(t.num_roots for t in tstemp.trees())
    recap_max_roots = max(t.num_roots for t in tsfull.trees())
    #print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
    #      f"After recapitation: {recap_max_roots}")

    numleaves = 0
    for k in firsttree.leaves():
        numleaves = numleaves + 1
    RANDOMSAMPLE = sample(range(numleaves), NUMBEROFLEAVESTOSAMPLE)
    #sample a random number of leaves to save.
    #CHECK THAT THERE IS ONLY 1 MUTATION
    numberofmutations = 0
    for k in tsfull.mutations():
        numberofmutations = numberofmutations + 1
        MutationNodeLocation = k.node
    if numberofmutations != 1:
        print("ERROR: NOT 1 MUTATION")

    IsDerived = [0] * len(RANDOMSAMPLE)
    for i in range(len(RANDOMSAMPLE)):
        currentnode = RANDOMSAMPLE[i]
        while currentnode >= 0:
            if currentnode == MutationNodeLocation:
                IsDerived[i] = 1
                break
            currentnode = firsttree.parent(currentnode)

    AncTimes = []
    DerTimes = []
    CROSSNODE = []
    for i in range(len(RANDOMSAMPLE)-1):
        for j in range(i + 1,len(RANDOMSAMPLE)):
            NodeTime = firsttree.mrca(RANDOMSAMPLE[i], RANDOMSAMPLE[j])
            if IsDerived[i] == 1 and IsDerived[j] == 1:
                DerTimes.append(NodeTime)
            elif IsDerived[i] == 0 and IsDerived[j] == 0:
                AncTimes.append(NodeTime)
            else:
                CROSSNODE.append(NodeTime)
                
    AncTimes = numpy.unique(AncTimes)
    DerTimes = numpy.unique(DerTimes)
    CROSSNODE = numpy.unique(CROSSNODE)
    #convert from nodes to times
    for i in range(len(AncTimes)):
        if AncTimes[i] < 0:
            print("PROBLEM2")
        else:
            AncTimes[i] = firsttree.time(AncTimes[i])
    for i in range(len(DerTimes)):
        DerTimes[i] = firsttree.time(DerTimes[i])
    for i in range(len(CROSSNODE)):
        CROSSNODE[i] = firsttree.time(CROSSNODE[i])
    AncTimes.sort()
    DerTimes.sort()
    CROSSNODE.sort()
    print(CROSSNODE[0])
    if len(DerTimes) + len(AncTimes) + 2 != NUMBEROFLEAVESTOSAMPLE:
        print("PROBLEM7")
    AncString = ""
    DerString = ""
    for time in AncTimes:
        AncString = AncString + str(time)  + ","
    for time in DerTimes:
        DerString = DerString + str(time) + ","
    f = open("INPUTTIMES/Times" + str(indexx) + ".txt", "w")
    indexx = indexx + 1
    f.writelines(DerString[0:(len(DerString)-1)]+ "\n" + AncString[0:(len(AncString)-1)] + "\n") # added indexing to remove comma at the end
    f.close()