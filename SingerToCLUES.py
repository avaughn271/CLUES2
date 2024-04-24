import numpy, argparse, tskit, os

"""Define the Arguments"""
parser = argparse.ArgumentParser()
parser.add_argument('--position',type=int, default=None)
parser.add_argument('--tree_path',type=str,default=None)
parser.add_argument('--output',type=str,default=None)

args =  parser.parse_args()
TREEPATH = args.tree_path

if os.path.exists(args.output + "_times.txt"):
  os.remove(args.output + "_times.txt")
STRINGLIST = os.listdir(TREEPATH)
boolean = [False] * len(STRINGLIST)

for i in range(len(STRINGLIST)):
    if len(STRINGLIST[i]) < 6:
        boolean[i] = False
    else:
        l = len(STRINGLIST[i])
        boolean[i] = ((STRINGLIST[i][(l - 6):]) == ".trees")
STRINGLIST = (numpy.array(STRINGLIST)[boolean]).tolist()

for k in STRINGLIST:
    POSITION = args.position

    mts = tskit.load(TREEPATH + "/" + k)

    firsttree = mts.first()
    numleaves = 0
    for k in firsttree.leaves():
        numleaves = numleaves + 1

    numtrees = 0
    for i in mts.trees():
        numtrees = numtrees + 1

    for i in range(len(mts.tables.sites)):
        if mts.tables.sites[i].position == POSITION:
            SITE = i
            break

    for i in range(len(mts.tables.mutations)):
        try:
            if mts.tables.mutations[i].site == SITE:
                NODE = mts.tables.mutations[i].node  # node above which mutation occurred.
                break
        except NameError:
            print("No mutation found at that position")
            exit()

    for i in mts.trees():
        if (i.interval.left <= POSITION and i.interval.right > POSITION):
            LOCALTREE = i
            break

    TotalNodes = []
    NodeArray = []
    for nod in LOCALTREE.nodes():
        TotalNodes.append(nod)
        NodeArray.append([nod, LOCALTREE.parent(nod)])
    NodeArray = numpy.array(NodeArray)

    def returnchildnode(currnode, NodeArray):
        if currnode not in NodeArray[:,1]:
            return [currnode]
        children = NodeArray[NodeArray[:,1] == currnode, 0]
        child1 = children[0]
        child2 = children[1]
        return [currnode] + returnchildnode(child1, NodeArray) + returnchildnode(child2, NodeArray)

    DerivedNodes = (returnchildnode(NODE, NodeArray) + [(NodeArray[NodeArray[:,0] == NODE, 1][0])] )
    AncestralNodes = numpy.setdiff1d(TotalNodes,DerivedNodes)

    DerivedTimes = []
    AncestralTimes = []
    for i in DerivedNodes:
        a = LOCALTREE.time(i)
        if a > 0:
            DerivedTimes.append(LOCALTREE.time(i))
    for i in AncestralNodes:
        a = LOCALTREE.time(i)
        if a > 0:
            AncestralTimes.append(LOCALTREE.time(i))
    DerivedTimes.sort()
    AncestralTimes.sort()

    AncString = ""
    DerString = ""
    for time in AncestralTimes:
        AncString = AncString + str(time)  + ","
    for time in DerivedTimes:
        DerString = DerString + str(time) + ","
    f = open(args.output + "_times.txt", "a")
    f.writelines(DerString[0:(len(DerString)-1)]+ "\n" + AncString[0:(len(AncString)-1)] + "\n") # added indexing to remove comma at the end
    f.close()