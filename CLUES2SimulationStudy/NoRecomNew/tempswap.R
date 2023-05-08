args = commandArgs(trailingOnly=TRUE)
L = as.numeric(args[1])
mu = as.numeric(args[2])
library(igraph)
iter = as.numeric(args[3])
thinningamount = as.numeric(args[4])
N = as.numeric(args[5]) #haploid popsize
output = toString(args[6])
inputnumber = toString(args[7])

A = read.table(paste0("INPUT/Data" , inputnumber ,".txt"))
mutations = A[-1,2]
mutationpatterns = table(mutations)
mutationcounts = unname(mutationpatterns)
mutationtypes = names(mutationpatterns)

Topology = data.matrix(read.table(paste0("INPUT/TreeTopology" , inputnumber ,".txt")) )
removenodes = as.numeric(names(table(Topology[,2] )[table(Topology[,2])== 1]))
print("removenodes")
print(removenodes)
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

###############################################################
Topology = Topology + 1 # this is now 1-indexed
if (nrow(Topology) != numberofbranches) {print("PROBLEM")}

Nodes = sort(unique(Topology[,1])) # no root node

TopologyConstraints = Topology
for (i in 1:(numberofleaves - 1)) {
  TopologyConstraints = rbind(TopologyConstraints, c(i,i+1))
}
for (i in 1:(numberofleaves - 1) ) {
  TopologyConstraints = rbind(TopologyConstraints, c(numberofleaves,i+numberofleaves))
}

NumberOfCounts = rep(0, length(Nodes))
ListOfLeafDescendants = list()
for (i in 1:(length(Nodes)+1) ) {
  ListOfLeafDescendants[[i]] = c(-1)
}
for (i in 1:numberofleaves) {
  currnode = i
  ListOfLeafDescendants[[currnode]] = c(ListOfLeafDescendants[[currnode]], i)
  while (currnode %in% Topology[,1] ) {
    currnode = Topology[which(Topology[,1] == currnode),2]

    ListOfLeafDescendants[[currnode]] = c(ListOfLeafDescendants[[currnode]],i)
  }
}

for (i in 1:length(ListOfLeafDescendants)) {
  ListOfLeafDescendants[[i]] = ListOfLeafDescendants[[i]][which(ListOfLeafDescendants[[i]] > -1)]
}

if (length(mutations) != 0) {

for (mutindex in 1:length(mutationcounts)) {
  currentmutationtype = mutationtypes[mutindex]
  numberofderived = lengths(gregexpr("1", currentmutationtype))
  firstmatch = gregexpr("1", currentmutationtype)[[1]]
  for (j in 1:length(ListOfLeafDescendants)) {
    if (length(ListOfLeafDescendants[[j]]) == length(firstmatch)) {
      if (sum(abs(firstmatch - ListOfLeafDescendants[[j]] )) == 0){
        NumberOfCounts[j] = mutationcounts[mutindex]
        break
      }
    }
  }
}
}
###############################  Permutations
checkpermutationoneorder = function(perm, order) { # order is vector of size 2.
  if (order[1] %in% perm & order[2] %in% perm ) {
    return(which(perm == order[1]) <= which(perm == order[2]) )
  }
  return(T)
}

checkpermutationallorders = function(perm, ordermatri) {
  for (i in 1:nrow(ordermatri)) {
    if (checkpermutationoneorder(perm, ordermatri[i,]) == F) {return(F)}
  }
  return(T)
}

getparent = function(top, child) {
  return(top[which(top[,1] == child),2]) #1 - indexed
}

###############################
print("generating random ordering")
#while (T) {
#  K = c(1:numberofleaves ,sample(
#    (numberofleaves + 1):(numberofleaves + numberofleaves-1),
#    numberofleaves - 1, replace = F   ))
#  if (checkpermutationallorders(K, TopologyConstraints)) break
#}

K = topo_sort(graph_from_data_frame(data.frame(TopologyConstraints)), mode = c("out"))

print(K)
print("done")
#checkorderings
CoalOrdering = order(K)
CoalOrdering = CoalOrdering[-length(CoalOrdering)]
NumberOfCounts = c(NumberOfCounts, 0)
Nodes = c(Nodes, length(NumberOfCounts))	
CoalOrdering = c(CoalOrdering, length(NumberOfCounts))
InterCoalTimes_notmapped = c(rep(0, numberofleaves), 100 * (1:(numberofleaves - 1)))
NodeTimes = rep(0, length(InterCoalTimes_notmapped))	 
branchlengths = rep(0, length(InterCoalTimes_notmapped))	 


CurrentTree = data.frame(NumberOfCounts, Nodes, CoalOrdering, InterCoalTimes_notmapped, NodeTimes, branchlengths)

for (i in (numberofleaves + 1):nrow(CurrentTree)) {
  CurrentTree$NodeTimes[which(CoalOrdering == i)] = (sum(CurrentTree$InterCoalTimes_notmapped[(numberofleaves + 1):i]))
}

for (i in 1:(nrow(CurrentTree) - 1) ) {
  CurrentTree$branchlengths[i] = CurrentTree$NodeTimes[Topology[which(Topology[,1] == i),2]] - CurrentTree$NodeTimes[i] 
}

logprior = function(tree) {
  coals = tree$InterCoalTimes_notmapped[(numberofleaves + 1):nrow(tree)]
  
  return(-sum(( choose(  rev(1 + (1:length(coals)  ) ) ,2  )  / N) * coals))
  
}

loglikelihood = function(tree) {
  total = 0
  branches = tree$branchlengths
  mutatincounts = tree$NumberOfCounts
  
  for (i in 1:(length(branches) - 1)) {
    total = total + dpois(mutatincounts[i], branches[i] * mu * L, log = TRUE)
  }
  return(total)
}

logposterior = function(tree) {return(logprior(tree) + loglikelihood(tree)   )}

adjustcoaltimes = function(tree) {
  indextopropose = sample((numberofleaves + 1):nrow(tree), 1)
  currentime = tree$InterCoalTimes_notmapped[indextopropose]
  tree$InterCoalTimes_notmapped[indextopropose] = rgamma(1, shape = 1,  scale = currentime) # rexp(1, 1/currentime) #changed
  newtime = tree$InterCoalTimes_notmapped[indextopropose]

  for (i in (numberofleaves + 1):nrow(tree)) {
    tree$NodeTimes[which(tree$CoalOrdering == i)] = (sum(tree$InterCoalTimes_notmapped[(numberofleaves + 1):i]))
  }
  
  for (i in 1:(nrow(tree) - 1) ) {
    tree$branchlengths[i] = tree$NodeTimes[Topology[which(Topology[,1] == i),2]] - tree$NodeTimes[i] 
  }

  return(list(tree, newtime, currentime))
}

adjustorderings = function(tree) {
  Rightbound = sample((numberofleaves + 1):(nrow(CurrentTree) - 1) ,1)
  Leftbound = Rightbound - 1
  rightindex = which(tree$CoalOrdering == Rightbound)
  leftindex = which(tree$CoalOrdering == Leftbound)
  temporary = tree$CoalOrdering
  temporary[leftindex] = Rightbound
  temporary[rightindex] = Leftbound
  if (!checkpermutationallorders(order(temporary), TopologyConstraints)) {
    return(list(tree, FALSE))
  }
  
  tree$CoalOrdering = temporary
  
  for (i in (numberofleaves + 1):nrow(CurrentTree)) {
    tree$NodeTimes[which(tree$CoalOrdering == i)] = (sum(tree$InterCoalTimes_notmapped[(numberofleaves + 1):i]))
  }
  
  for (i in 1:(nrow(tree) - 1) ) {
    tree$branchlengths[i] = tree$NodeTimes[Topology[which(Topology[,1] == i),2]] - tree$NodeTimes[i]
  }
  return(list(tree, TRUE))
}

acc = rep(0,length((numberofleaves + 1):nrow(CurrentTree)))
swaps = 0
scaleattempts = rep(0,length((numberofleaves + 1):nrow(CurrentTree)))
swapattempts = 0
thinningsets = round(seq(1,iter, length.out = thinningamount + 1))[-1]

print("Start")


for (i in 1:iter) {

if (i %% 2000 == 0) {print(i)}

  if (runif(1) < 0.5) {
    proposedtree = data.frame(CurrentTree)
    unpack = adjustcoaltimes(proposedtree)
    proposedtree = unpack[[1]]
    newtime = unpack[[2]]
    oldtime = unpack[[3]]

    logmhratio  = dgamma(oldtime, shape = 1, scale = newtime, log = T) - dgamma(newtime, shape = 1, scale = oldtime, log = T) #log(oldtime/newtime) - oldtime/newtime + newtime/oldtime
    changedindex = (which(CurrentTree[[4]] != proposedtree[[4]]) - numberofleaves)
    if (log(runif(1)) <=  logmhratio + logposterior(proposedtree) - logposterior(CurrentTree)  ) {
      CurrentTree = proposedtree
      acc[changedindex] = acc[changedindex]  + 1
    }
    scaleattempts[changedindex] =  scaleattempts[changedindex] + 1
  }
  else {

    swapattempts = swapattempts + 1
    proposedtree = data.frame(CurrentTree)
    unpack = adjustorderings(proposedtree)

    proposedtree = unpack[[1]]
    if (unpack[[2]]){
    if (log(runif(1)) <=  loglikelihood(proposedtree) - loglikelihood(CurrentTree)  ) { # cahnged to only log likelihood is changing
      CurrentTree = proposedtree
      swaps = swaps  + 1
    }
    }

    }
  if (i %in% thinningsets) {
  write.csv(CurrentTree, paste0(output, "/results", toString(i), ".csv") , row.names=FALSE)

  }
}
print("Scale acceptances")
print(acc/scaleattempts)
print("Swap acceptances")

print(swaps/swapattempts)
