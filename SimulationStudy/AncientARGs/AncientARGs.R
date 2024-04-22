args = commandArgs(trailingOnly=TRUE)
N = as.numeric(args[1]) # this is haploid Ne
selectioncoef = as.numeric(args[2])
h = 0.5 ##1, 1 + hs, 1 + s
GenerationsBetweenSamples = args[3]
mincurrfrequnecy = as.numeric(args[4])
DerivedTimes = c()
AncestralTimes = c()

MAXTIME_BeforePresent =  as.numeric(args[5])
numsampledlineages = as.numeric(args[6])

OneTrajectoryFunction = function() {
  while(T) {
    timeofmutation = sample(1:MAXTIME_BeforePresent, 1)
    #print(timeofmutation)
    TotalFreq = rep(-1, timeofmutation)
    TotalFreq[1] = 1/N # This is the initial frequency
    currentfreq = TotalFreq[1]
    nextgeneration = 2
    while(currentfreq > 0 & currentfreq < 1) {
      p = currentfreq
      newmean = p - (selectioncoef* (-1 + p)* p* (-p + h *(-1 + 2* p))) /(-1 + selectioncoef *(2* h *(-1 + p) - p) *p)
      currentfreq = rbinom(1, N, newmean) / N
      TotalFreq[nextgeneration] = currentfreq
      nextgeneration = nextgeneration + 1
      if (  (nextgeneration > length(TotalFreq)) &    (currentfreq > mincurrfrequnecy) &   (currentfreq < 0.96)){ # enforce that modern freq is at least 25%.
        print(c(timeofmutation, currentfreq))
        return(TotalFreq)
      }
    }
  }
}
print("AAA")
OneTrajectory = rev(OneTrajectoryFunction())
OneTrajectory = c(OneTrajectory, rep(0,100000))
currderivednumber = round(numsampledlineages * OneTrajectory[1])
currancestralnumber = round(numsampledlineages - currderivednumber)

if (currderivednumber >= numsampledlineages - 2) {
  currderivednumber =  numsampledlineages - 3
  currancestralnumber = round(numsampledlineages - currderivednumber)
}

if (currancestralnumber >= numsampledlineages - 2) {
  currancestralnumber =  numsampledlineages - 3
  currderivednumber = round(numsampledlineages - currancestralnumber)
}

DerivedCoals = c()
AncestralCoals = c()

for (ii in 1:(10000 + 3) ) {
  currentfreq = OneTrajectory[ii+1]

  if (round(currentfreq * N) == 0 & currderivednumber > 0) {
    currderivednumber = currderivednumber - 1
    currancestralnumber = currancestralnumber + 1
  }
  if (currderivednumber > 0) {
parents = sample(round(currentfreq * N),  currderivednumber  , replace = T)
coals = length(parents) - length(unique(parents))
DerivedCoals = c(DerivedCoals, rep(ii, coals))
currderivednumber = currderivednumber - coals
}

parents = sample(round(N - currentfreq * N),  currancestralnumber  , replace = T)
coals = length(parents) - length(unique(parents))
AncestralCoals = c(AncestralCoals, rep(ii, coals))
currancestralnumber = currancestralnumber - coals
#print(c(length(parents) , length(unique(parents)), coals))

if (ii == 100 | ii == 50) {
  derr = sample(0:1, size = 160, prob = c(1 - currentfreq, currentfreq), replace = T)
nderrr = sum(derr)
nanccc = length(derr) - sum(nderrr)
currderivednumber = currderivednumber + nderrr
currancestralnumber = currancestralnumber + nanccc
DerivedTimes = c(DerivedTimes, rep(ii, nderrr))
AncestralTimes = c(AncestralTimes, rep(ii, nanccc))
}

}

#print(DerivedTimes)
#print(AncestralTimes)
#print(AncestralCoals)
#print(DerivedCoals)

dertimess =  toString(DerivedTimes[1])
for (i in 2:length(DerivedTimes)) {
  dertimess = paste0(dertimess, ";", toString(DerivedTimes[i]))
}

if (length(AncestralTimes) == 0 ) {
  anctimess = ""
} else {
anctimess =  toString(AncestralTimes[1])
for (i in 2:length(AncestralTimes)) {
  anctimess = paste0(anctimess, ";", toString(AncestralTimes[i]))
}
}

writeLines(toString(OneTrajectory[1]), "ModernFreq.txt", sep = "\n")
currderivednumber = round(numsampledlineages * OneTrajectory[1])
currancestralnumber = round(numsampledlineages - currderivednumber)

if (currderivednumber >= numsampledlineages - 2) {
  currderivednumber =  numsampledlineages - 3
  currancestralnumber = round(numsampledlineages - currderivednumber)
}

if (currancestralnumber >= numsampledlineages - 2) {
  currancestralnumber =  numsampledlineages - 3
  currderivednumber = round(numsampledlineages - currancestralnumber)
}

print(c(     currderivednumber ,  length(DerivedCoals)    , length(DerivedTimes)            ))
DerivedCoals = c(DerivedCoals, rep(200000, currderivednumber- length(DerivedCoals) + length(DerivedTimes)))
AncestralCoals = c(AncestralCoals, rep(200000, currancestralnumber - 1 - length(AncestralCoals)  + length(AncestralTimes) ))

derstring = toString(DerivedCoals[1])
for (i in 2:length(DerivedCoals)) {
  derstring = paste0(derstring, ",", toString(DerivedCoals[i]))
}
ancstring = toString(AncestralCoals[1])
for (i in 2:length(AncestralCoals)) {
  ancstring = paste0(ancstring, ",", toString(AncestralCoals[i]))
}
writeLines( c(dertimess, anctimess,derstring , ancstring),  "INPUTTIMES/AncientTimes.txt")
