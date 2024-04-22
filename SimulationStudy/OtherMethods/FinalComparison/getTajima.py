import msprime, sys

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
    A = msprime.load("../../NoRecomNewer/plotG/val" + str(val) + "/INPUT/Tree" + str(i + 1) + ".trees")
    
    with open("TrueResults/Tajima" + str(i) +".txt", 'w') as fp:
        fp.write("%s\n" % A.Tajimas_D())