args = commandArgs(trailingOnly=TRUE)

get01fromtable = function(tablee, topology, leaf1, leaf2) {
  parents0 = c()
  currnode = leaf1
  while (currnode != 0) {
    parents0 = c(parents0 , currnode)
    currnode = topology[[2]][which(topology[[1]] == currnode)]
  }
  parents1 = c()
  currnode = leaf2
  while (currnode != 0) {
    parents1 = c(parents1 , currnode)
    currnode = topology[[2]][which(topology[[1]] == currnode)]
  }
  commonparents = intersect(parents1,parents0)
  
  rr=( min(tablee$NodeTimes[min(commonparents)]) )
  return(rr)
}

RelevantFolder = paste0("INPUT/samples", toString(args[1]))

IsDerived = readLines(paste0("INPUT/IsDerived", toString(args[1]), ".txt"))

ListOfFiles = list.files(RelevantFolder, pattern = "*.csv")
Samples = c()
Topology = read.table(paste0("Input/TreeTopology", toString(args[1]), ".txt"), sep ="\t") + 1

#####

Topology = data.matrix(read.table(paste0("INPUT/TreeTopology", toString(args[1]) ,".txt")) )
removenodes = as.numeric(names(table(Topology[,2] )[table(Topology[,2])== 1]))
removenodes = removenodes[removenodes > 0]
if (length(removenodes) > 0) {
  for (i in removenodes) {
    middlenumber = Topology[which(Topology[,2] == i),2]
    lowernumber = Topology[which(Topology[,2] == i),1]
    topnumber = Topology[which(Topology[,1] == middlenumber),2]
    Topology[which(Topology[,1] == lowernumber),2] = topnumber
    Topology = Topology[-which(Topology[,1] == middlenumber),]
  }
}
Topology = Topology[-which(Topology[,2]  < 0 ),]

numberofleaves = length(setdiff(Topology[,1], Topology[,2]))

numberofbranches = numberofleaves + numberofleaves - 2

desirednodes = (1:length(unique(c(Topology))))-1
actualnodes = unique(c(Topology))
desirednodesnew = setdiff(desirednodes,  intersect(desirednodes,actualnodes))
actualnodesnew = setdiff(actualnodes,  intersect(desirednodes,actualnodes))
for (i in 1:length(desirednodesnew)) {
  Topology[which(Topology == actualnodesnew[i])] = desirednodesnew[i]
}
Topology = data.frame(Topology + 1)
Topology = rbind(Topology, c(max(Topology[[2]]), 0  ))
####
for (k in ListOfFiles) {
TreeFile = read.csv(paste0(RelevantFolder, "/",k), sep =",")

DerivedCoals = c()
AncestralCoals = c()
for (leaf1 in 1:(length(IsDerived)-1) ) {
  for (leaf2 in (leaf1+1):length(IsDerived)) {
    if (IsDerived[leaf1] == "1" & IsDerived[leaf2] == "1") {
      DerivedCoals = c(DerivedCoals, get01fromtable(TreeFile, Topology, leaf1, leaf2) )
    } else if (IsDerived[leaf1] == "0" & IsDerived[leaf2] == "0") {
      AncestralCoals = c(AncestralCoals, get01fromtable(TreeFile,Topology, leaf1, leaf2) )
    }
  }
}
DerivedCoals = unique(DerivedCoals)
AncestralCoals = unique(AncestralCoals)

if (length(DerivedCoals) + length(AncestralCoals) != length(IsDerived) - 2) {
  print("INFINITE SITES ASSUMPTION NOT SATISFIED")
}
DerivedCoals = c(DerivedCoals, setdiff((TreeFile$NodeTimes)[TreeFile$NodeTimes > 0], c(AncestralCoals, DerivedCoals)))

if (length(DerivedCoals) + length(AncestralCoals) != length(IsDerived)  - 1) {
  print("INFINITE SITES ASSUMPTION NOT SATISFIED")
}

anc = sort(AncestralCoals)
der = sort(DerivedCoals)
write(c(paste(der, collapse = ","), paste(anc, collapse = ",")), paste0("INPUTTIMES/Times", toString(args[1]),".txt" ),append=TRUE)
}
