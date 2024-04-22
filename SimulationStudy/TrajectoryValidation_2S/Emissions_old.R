N = 30000

TCutoff = 400
TRIALS = 1000
h = 0.5
MAXTIME_BeforePresent = 10000
OneTrajectoryFunction = function(selectioncoef1, selectioncoef2) {
  while(T) {
    timeofmutation = sample(1:MAXTIME_BeforePresent, 1)
    TotalFreq = rep(-1, timeofmutation)
    TotalFreq[1] = 1/N # This is the initial frequency
    currentfreq = TotalFreq[1]
    nextgeneration = 2
    while(currentfreq > 0 & currentfreq < 1) {
      p = currentfreq
      if (nextgeneration > (timeofmutation -  200    ))  {selectioncoef = selectioncoef2} else 
        {selectioncoef = selectioncoef1}
      newmean = p - (selectioncoef* (-1 + p)* p* (-p + h *(-1 + 2* p))) /(-1 + selectioncoef *(2* h *(-1 + p) - p) *p)
      currentfreq = rbinom(1, N, newmean) / N
      TotalFreq[nextgeneration] = currentfreq
      nextgeneration = nextgeneration + 1
      if ( (nextgeneration > length(TotalFreq)) & (abs(TotalFreq[nextgeneration - 1] - 0.25 ) < 0.01)) {
     #   print(timeofmutation)
        return(TotalFreq)
    }
    }}
    print("STOP")
}

accvals = 0
while (T) {
  selcoef1 = runif(1, 0.01, 0.2)
  selcoef2 = runif(1, -0.08, 0)
  
 # print(c("trying with ", selcoef1,selcoef2))
SampledFrequency = c(rev(OneTrajectoryFunction(selcoef1,selcoef2)), rep(0,1000))
#print("done")
retainboolean = T
for (ii in 1:(TCutoff + 3) ) {
if ( ii %in% c(100,200,250,350) ) {
  currentfreq = SampledFrequency[ii]
  r = runif(1)
  if (r < currentfreq**2) {coll = 2}
  else if (r < currentfreq**2 + (1 - currentfreq)**2 ) { coll = 0}
  else {coll = 1}
  if (coll == 2 &  ii != 200 ) {retainboolean = F; break}
#  print("place1")
  
  if (coll == 0 & ii  != 350 )  {retainboolean = F; break}
#  print("place2")
  
  if (coll == 1 &   ii %in% c(200, 350)   )  {retainboolean = F; break}
#  print("place3")
  
}

}
if (retainboolean) {
  print(c(selcoef1,selcoef2))
  accvals = accvals + 1
  write.table(SampledFrequency,paste0("Frequencies/Freqs", toString(accvals) ,".txt"), row.names = F, col.names = F)
}
if (accvals == TRIALS) { break}
}