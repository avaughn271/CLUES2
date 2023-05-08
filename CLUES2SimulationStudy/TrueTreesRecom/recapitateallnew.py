import numpy, msprime, pyslim, sys

HaploidNe = float(sys.argv[1])
NUMBEROFLEAVESTOSAMPLE = int(sys.argv[2])
sss = float(sys.argv[3])
numtrials = int(sys.argv[4])
LL=1e7
if sss == 0.0: sss = 1e-6 # so that the code works properly.

indexx = 1
for filename in range(numtrials):
    # define hard sweep model
    sweep_model = msprime.SweepGenicSelection(
        position= LL/2,  # middle of chrom
        start_frequency=1.0 / (2 * HaploidNe),
        end_frequency=0.75,
        s=(sss*2),
        dt=1e-6)

    reps = msprime.sim_ancestry(
        NUMBEROFLEAVESTOSAMPLE/2,
        model=[sweep_model],
        population_size=HaploidNe/2.0,
        recombination_rate=1e-8,
        sequence_length=LL)

    demm = msprime.Demography.from_tree_sequence(reps)
    for pop in demm.populations:
        pop.initial_size = HaploidNe / 2.0  ##POSITON1 This is diploid size

    tsfull = pyslim.recapitate(reps, demography= demm, recombination_rate=0.0)

    firsttree = reps.first() #CHANGE!!!!!!!!!
    for t in reps.trees():
        if t.interval[0] <= LL/2 and t.interval[1] > LL/2:
            firsttree = t
            print("done1")
            break
    counttt = 0 
    for t in reps.trees():
        counttt = counttt + 1

    print("Number of trees is:" + str(counttt))

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
    firsttree = tsfull.first()  #CHANGE!!!!!!!!!
    for t in tsfull.trees():
        if t.interval[0] <= LL/2 and t.interval[1] > LL/2:
            firsttree = t
            print("done2")
            break

    AncTimes = []
    DerTimes = []
    CROSSNODE = []
    for i in range(len(IsDerived)-1):
        for j in range(i + 1,len(IsDerived)):
            NodeTime = firsttree.mrca(i, j)
            if IsDerived[i] == 1 and IsDerived[j] == 1:
                DerTimes.append(NodeTime)
            elif IsDerived[i] == 0 and IsDerived[j] == 0:
                AncTimes.append(NodeTime)
            else:
                CROSSNODE.append(NodeTime)
                
    AncTimes = numpy.unique(AncTimes)
    DerTimes = numpy.unique(DerTimes)
    CROSSNODE = numpy.unique(CROSSNODE)
    NewAncTimes = [0.0] * len(AncTimes)
    NewDerTimes = [0.0] * len(DerTimes)
    NewCROSSNODE = [0.0] * len(CROSSNODE)

    #convert from nodes to times
    for i in range(len(AncTimes)):
        if AncTimes[i] < 0:
            print("PROBLEM2")
        else:
            NewAncTimes[i] = float(firsttree.time(AncTimes[i]))
    for i in range(len(DerTimes)):
        NewDerTimes[i] = float(firsttree.time(DerTimes[i]))
    for i in range(len(CROSSNODE)):
        NewCROSSNODE[i] = float(firsttree.time(CROSSNODE[i]))
    NewDerTimes.append(NewCROSSNODE[0])
    NewAncTimes.sort()
    NewDerTimes.sort()
    if len(NewDerTimes) + len(NewAncTimes) + 1 != NUMBEROFLEAVESTOSAMPLE:
        print("PROBLEM7")
    AncString = ""
    DerString = ""
    for time in NewAncTimes:
        AncString = AncString + str(time)  + ","
    for time in NewDerTimes:
        DerString = DerString + str(time) + ","
    f = open("INPUTTIMES/Times" + str(indexx) + ".txt", "w")
    indexx = indexx + 1
    f.writelines(DerString[0:(len(DerString)-1)]+ "\n" + AncString[0:(len(AncString)-1)] + "\n") # added indexing to remove comma at the end
    f.close()