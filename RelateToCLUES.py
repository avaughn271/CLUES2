from Bio import Phylo
from io import StringIO
import pandas as pd
import re
import numpy as np
import argparse

def parse_args():
	"""Define the Arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--RelateSamples',type=str, default=None)
	parser.add_argument('--DerivedFile',type=str, default=None)
	parser.add_argument('--out',type=str,default=None)
	return parser.parse_args()

def newick_to_list(newick_string):
    name_regex = re.compile(r"[^(),:]+")  # Matches node names

    def traverse(node, parent=None):
        if node is None:
            return
        for child in node:
            if isinstance(child, str):
                df.loc[len(df)] = [parent, child]
            else:
                name_match = name_regex.search(child.name)
                node_name = name_match.group() if name_match else None
                df.loc[len(df)] = [parent, node_name]
                traverse(child.clades, parent=node_name)

    df = pd.DataFrame(columns=["parent", "child"])

    tree = Phylo.read(StringIO(newick_string), "newick")

    traverse(tree.clade, parent=None)

    Parents = []
    Childs = []
    p = (df.iloc[:,0])
    c = (df.iloc[:,1])
    for i in p: Parents.append(i)
    for i in c: Childs.append(i)
    for k in range(len(Parents)):
        if Parents[k] == None:
            Parents[k] = "root"
    return([Childs, Parents])

def CalculateDerivedBelow(List, nodee, nodeeindex):
    if not (nodee in List[1]):
        List[3][nodeeindex] = List[2][nodeeindex]
    else:
        child1index = (np.where(np.array(List[1]) == nodee))[0][0]
        child2index = (np.where(np.array(List[1]) == nodee))[0][1]
        CalculateDerivedBelow(List, List[0][child1index], child1index)
        CalculateDerivedBelow(List, List[0][child2index], child2index)
        List[3][nodeeindex] = List[3][child1index] + List[3][child2index]

def CalculateAncestralBelow(List, nodee, nodeeindex):
    if not (nodee in List[1]):
        List[4][nodeeindex] = 1 - List[2][nodeeindex]
    else:
        child1index = (np.where(np.array(List[1]) == nodee))[0][0]
        child2index = (np.where(np.array(List[1]) == nodee))[0][1]
        CalculateAncestralBelow(List, List[0][child1index], child1index)
        CalculateAncestralBelow(List, List[0][child2index], child2index)
        List[4][nodeeindex] = List[4][child1index] + List[4][child2index]

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
    IsDerived = isder
    f = open(args.RelateSamples, "r")

    TotalStrings = []
    Lines = (f.readlines())[1:]
    newicktree = (Lines[0].split("\t"))
    newicktree = newicktree[len(newicktree)-1]
    newicktree = re.sub(r':\d+(\.\d+)?,', ',', newicktree)
    newicktree = re.sub(r':\d+(\.\d+)?\)', ')', newicktree)

    letter = "a"
    alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    def nextletter(let):
        lastletter = let[len(let) - 1]
        if lastletter != "z":
            ii = alphabet.index(lastletter)
            return(let[0:(len(let)-1)] + alphabet[ii+1])
        else:
            return(let + "a")

    while True:
        try:
            aaaa = (newicktree.index("),"))
            newicktree = newicktree[0:aaaa]  + ")" + letter + ","  + newicktree[(aaaa+2):len(newicktree)]
            letter = nextletter(letter)
        except ValueError:
            try:
                bbbb = (newicktree.index("))"))
                newicktree = newicktree[0:bbbb] + ")" + letter + ")" + newicktree[(bbbb+2):len(newicktree)]
                letter = nextletter(letter)
            except ValueError:
                break

    List = newick_to_list(newicktree)

    TotalDerived = np.sum(IsDerived)
    TotalAncestral = len(IsDerived) - TotalDerived
    #remove branch length by substituting :____, by , nad replacing :____) by )
    #Then, add a letter after each closed parenthesis.
    List.append([-1] * len(List[0]))

    nnnnn = np.array(List[0])
    for i in range(len(IsDerived)):
        List[2][(np.where(nnnnn == str(i)))[0][0] ] = IsDerived[i]

    List.append ( [0] * len(List[0])) # number of derived below
    List.append ( [0] * len(List[0])) # number of ancestral below

    child1index = (np.where(np.array(List[1]) == "root"))[0][0]
    child2index = (np.where(np.array(List[1]) == "root"))[0][1]
    CalculateDerivedBelow(List,  List[0][child1index], child1index)
    CalculateDerivedBelow(List,  List[0][child2index], child2index)

    child1index = (np.where(np.array(List[1]) == "root"))[0][0]
    child2index = (np.where(np.array(List[1]) == "root"))[0][1]
    CalculateAncestralBelow(List,  List[0][child1index], child1index)
    CalculateAncestralBelow(List,  List[0][child2index], child2index)

    TotalError = [0] * len(List[0])
    for i in range(len(TotalError)):
        TotalError[i] = TotalDerived - List[3][i] + List[4][i]
    minflips = min(TotalError)
    if minflips == 0:
        print("Infinite sites assumption satisfied. No allele flips necessary.")
    else:
        minnn = min(TotalError)
        if minnn == 1:
            print("Infinite sites assumption not satisfied. Flipping " + str(int(minnn)) + " leaf out of "  + str(len(IsDerived)) + " leaves.")
        else:
            print("Infinite sites assumption not satisfied. Flipping " + str(int(minnn)) + " leaves out of " + str(len(IsDerived)) + " leaves.")
        FlippingIndex = TotalError.index(minnn)

        NewIsDerived = [0] * len(IsDerived)
        Descendants = []

        def FindDescendants(List, nodee, nodeeindex):
            if not (nodee in List[1]):
                Descendants.append(int(List[0][nodeeindex]))
            else:
                child1index = (np.where(np.array(List[1]) == nodee))[0][0]
                child2index = (np.where(np.array(List[1]) == nodee))[0][1]
                FindDescendants(List, List[0][child1index], child1index)
                FindDescendants(List, List[0][child2index], child2index)

        FindDescendants(List,  List[0][FlippingIndex], FlippingIndex)

        for i in Descendants:
            NewIsDerived[i] = 1
        for i in range(len(NewIsDerived)):
            if NewIsDerived[i] > IsDerived[i]:
                print("Flipped leaf " + str(i) + " from ancestral (0) to derived (1).")
            if NewIsDerived[i] < IsDerived[i]:
                print("Flipped leaf " + str(i) + " from derived (1) to ancestral (0).")
        
    for i in range(len(Lines)):
        newicktree = (Lines[i].split("\t"))
        newicktree = newicktree[len(newicktree)-1]
        TotalStrings.extend(OneTreeToList2((newicktree[:-1]), NewIsDerived))
    
    if minflips > 0:
        tempderived = []
        for i in NewIsDerived: tempderived.append(str(i) + "\n")
        f = open(args.out + "_derived.txt", "w+")
        f.writelines(tempderived)
        f.close()
    f = open(args.out + "_times.txt", "w+")
    f.writelines(TotalStrings)
    f.close()