import msprime, sys

def getH12(treeseq):
    Haplotypes = [[0]]
    bool = True
    for i in range(1,int(float(sys.argv[2]))):
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
    with open("TrueResults/H12_" + str(i) +".txt", 'w') as fp:
        fp.write("%s\n" % getH12(A))