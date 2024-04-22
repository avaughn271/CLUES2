import msprime, sys
import numpy as np
LLL = 1e6
nummleavess = int(float(sys.argv[2]))

def getnSL(treeseq, isderived):
    numhaps = len(isderived)
    H = np.zeros((numhaps, numhaps))
    GenotypeMatrix = treeseq.genotype_matrix()
    Positions = []
    for i in range(len(treeseq.tables.sites)):
        Positions.append(treeseq.tables.sites[i].position)
    for i in range(1,len(Positions)):
        if Positions[i - 1] <= LLL/2 and Positions[i] >= LLL/2:
            origleftt = i - 1
            origrightt = i
            break
    for i in range(numhaps):
        for j in range(numhaps):
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
    nder = np.sum(isderived)
    nanc = len(isderived) - np.sum(isderived)
    SLA = 0
    SLD = 0
    for i in range(numhaps):
        for j in range(numhaps):
            if i < j:
                if isderived[i] == 1 and isderived[j] == 1:
                    SLD = SLD + H[i,j]
                if isderived[i] == 0 and isderived[j] == 0:
                    SLA = SLA + H[i,j]
    SLD = SLD * 2 / (nder * (nder - 1))
    SLA = SLA * 2 / (nanc * (nanc - 1))
    return(np.log(SLA/SLD))

sss = float(sys.argv[1])
if sss < 0.00001:
    val = 1
elif sss < 0.00101:
    val = 2
elif sss < 0.0025001:
    val = 3
elif sss < 0.005001:
    val = 4
elif sss < 0.00750001:
    val = 5
elif sss < 0.01001:
    val = 6

for i in range(30):
    A = msprime.load("../INPUT" + str(val) + "/Tree" + str(i + 1) + ".trees")
    my_file = open("../INPUT" + str(val) + "/IsDerived" + str(i + 1) + ".txt")
    
    data = my_file.read()
    
    IsDerived = [eval(ivvvv) for ivvvv in (data.split("\n"))[0:24]]
    my_file.close()
    with open("TrueResults/nSL_" + str(i) +".txt", 'w') as fp:
        fp.write("%s\n" % getnSL(A, IsDerived))