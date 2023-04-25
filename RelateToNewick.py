from Bio import Phylo
from io import StringIO
import numpy as np
import argparse

def parse_args():
	"""Define the Arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--RelateSamples',type=str, default=None)
	parser.add_argument('--DerivedFile',type=str, default=None)
	parser.add_argument('--out',type=str,default=None)
	return parser.parse_args()

def NewickToTMRCAMAtrix(newickstring):
    # First, read in the tree
    tree = Phylo.read(StringIO(newickstring), 'newick')
    # Get the leaves
    leaves = tree.get_terminals()
    L = len(leaves)
    # Create a L by L matrix to store the pairwise TMRCA
    TMRCA = np.zeros((L, L))
    # For each pair of leaves, get the TMRCA
    MaxDepth = -1
    for i in tree.depths().keys():
        if tree.depths()[i] > MaxDepth:
            MaxDepth = tree.depths()[i] 

    LeafNames = []
    for i in range(L):
        LeafNames.append(int(leaves[i].name))
    for i in range(L-1):
        for j in range(i+1, L):
            TMRCA[i, j] = MaxDepth - ((tree.depths())[tree.common_ancestor(leaves[i].name, leaves[j].name)])
            #TMRCA[i, j] = round(tree.distance(leaves[i], leaves[j]) / 2.0 , 3) # round to account for rounding problems
    # Return the matrix
    return TMRCA, LeafNames

def converttofile(TMRCA ,LeafNames, IsDerived):
    numleafs = len(LeafNames)
    Derivedbranches = []
    Ancestralbrnaches = []
    Mixed = []
    for leaf1 in range(numleafs-1):
        for leaf2 in range(leaf1 + 1, numleafs):
            if IsDerived[LeafNames[leaf1]] == 1 and IsDerived[LeafNames[leaf2]] == 1:
                Derivedbranches.append(TMRCA[leaf1, leaf2])
            elif IsDerived[LeafNames[leaf1]] == 0 and IsDerived[LeafNames[leaf2]] == 0:
                Ancestralbrnaches.append(TMRCA[leaf1, leaf2])
            else:
                Mixed.append(TMRCA[leaf1, leaf2])

    Ancestralbrnaches = (np.unique(Ancestralbrnaches))
    Derivedbranches = (np.unique(Derivedbranches))
    Mixedtemp = (np.unique(Mixed))
    if not (len(Ancestralbrnaches) + len(Derivedbranches) == numleafs - 2):
        print("Infinite sites assumption not satisfied. Incorrect output will follow.")
    Mixed = []
    for element in Mixedtemp:
        if element not in Ancestralbrnaches and element not in Derivedbranches:
            Mixed.append(element)
    if not (len(Ancestralbrnaches) + len(Derivedbranches) + len(Mixed) == numleafs - 1):
        print("Infinite sites assumption not satisfied. Incorrect output will follow.")
    return(np.sort(Ancestralbrnaches), np.sort(np.concatenate((Derivedbranches, Mixed))))

#AA =  '(5:6.0,(4:5.0,(3:4.2,(1:2.5,(2:2,0:2.000000001):0.5):1.7):0.8):1.0);' ###assumes input is 0-indexed
#BB = [1,1,1,0,0,0]

def OneTreeToList2(Newick,IsDerived):
    A = NewickToTMRCAMAtrix(Newick)
    Ancestralbrnaches,Derivedbranches = converttofile(A[0], A[1], IsDerived)
    AncString = ""
    DerString = ""
    for time in Ancestralbrnaches:
        AncString = AncString + str(time)  + ","
    for time in Derivedbranches:
        DerString = DerString + str(time) + ","
    DerString=DerString[:-1] + "\n"
    AncString=AncString[:-1] + "\n"
    return(DerString, AncString)

if __name__ == "__main__":
    args = parse_args()                    
    isder = np.loadtxt(args.DerivedFile)
    f = open(args.RelateSamples, "r")

    TotalStrings = []
    Lines = (f.readlines())[1:]
    for i in range(len(Lines)):
        newicktree = (Lines[i].split("\t"))
        newicktree = newicktree[len(newicktree)-1]
        TotalStrings.extend(OneTreeToList2((newicktree[:-1]), isder))

    f = open(args.out + ".txt", "w+")
    f.writelines(TotalStrings)
    f.close()
