import numpy, tskit
from random import sample

tsfull = tskit.load("slim.trees")
a = tsfull.first()
NUMBEROFLEAVESTOSAMPLE = 100

numleaves = 0
for k in a.leaves():
    numleaves = numleaves + 1
RANDOMSAMPLE = sample(range(numleaves),   NUMBEROFLEAVESTOSAMPLE)

#CHECK THAT THERE IS ONLY 1 TREE AND ONLY 1 MUTATION and 1 single root
numberoftrees = 0
for numtrees in tsfull.trees():
    numberoftrees = numberoftrees + 1
if numberoftrees != 1:
    print("ERROR: WRONG NUMBER OF TREES")
numberofmutations = 0
for k in tsfull.mutations():
    numberofmutations = numberofmutations + 1
    MutationNodeLocation = k.node
if numberofmutations != 1:
    print("ERROR: NOT 1 MUTATION")
if not a.has_single_root:
    print("ERROR: NOT 1 ROOT")
###################################################

IsDerived = [0] * len(RANDOMSAMPLE)
for i in range(len(RANDOMSAMPLE)):
    currentnode = RANDOMSAMPLE[i]
    while currentnode >= 0:
        if currentnode == MutationNodeLocation:
            IsDerived[i] = 1
            break
        currentnode = a.parent(currentnode)

AncTimes = []
DerTimes = []

for i in range(len(RANDOMSAMPLE)-1):
    for j in range(i + 1,len(RANDOMSAMPLE)):
        NodeTime = a.mrca(RANDOMSAMPLE[i], RANDOMSAMPLE[j])
        if IsDerived[i] == 1 and IsDerived[j] == 1:
            DerTimes.append(NodeTime)
        else:
            if IsDerived[i] + IsDerived[j] == 1 and NodeTime > 0:
                print("MIXES")
            AncTimes.append(NodeTime)
AncTimes = numpy.unique(AncTimes)
DerTimes = numpy.unique(DerTimes)

#convert from nodes to times
for i in range(len(AncTimes)):
    if AncTimes[i] < 0:
        AncTimes[i] = 10000000.0
    else:
        AncTimes[i] = a.time(AncTimes[i])
for i in range(len(DerTimes)):
    DerTimes[i] = a.time(DerTimes[i])

AncTimes.sort()
DerTimes.sort()
AncTimes = AncTimes[0:(len(AncTimes) - 1)]
if len(AncTimes) == 0 or len(DerTimes) == 0:
    print("ERROR: NO VARIATION IN GENOTYPES")
for i in range(NUMBEROFLEAVESTOSAMPLE + 3 - len(AncTimes) - len(DerTimes)):
    AncTimes = numpy.append(AncTimes, 1000000.0)

AncString = ""
DerString = ""
print(DerTimes)
A = (numpy.sum(IsDerived) - 1)
B = (NUMBEROFLEAVESTOSAMPLE - numpy.sum(IsDerived) - 1)
for time in AncTimes[0:A]:
    AncString = AncString + str(time)  + ","
for time in DerTimes[0:B]:
    DerString = DerString + str(time) + ","
f = open("Times.txt", "w")

f.writelines(DerString[0:(len(DerString)-1)]+ "\n" + AncString[0:(len(AncString)-1)] + "\n") # added indexing to remove comma at the end
f.close()