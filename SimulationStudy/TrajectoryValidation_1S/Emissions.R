N = 30000

TCutoff = 400
TRIALS = 50
h = 0.5
MAXTIME_BeforePresent = 10000
OneTrajectoryFunction = function(selectioncoef) {
  while(T) {
    timeofmutation = sample(1:MAXTIME_BeforePresent, 1)
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
      if ( nextgeneration > length(TotalFreq) ) {
        print(timeofmutation)
        return(TotalFreq)
    }
    }}
    return(c(-1,-1))
}

accvals = 0
while (T) {
retainboolean = T
  selcoef = runif(1, 0.001, 0.1)
  print(c("trying with ", selcoef))
SampledFrequency = c(rev(OneTrajectoryFunction(selcoef))     , rep(0,1000))
print(SampledFrequency[1])
if (abs(SampledFrequency[1] - 0.75 ) > 0.01) {next}
print("done")
for (ii in 1:(TCutoff + 3) ) {
if ( (ii  %% 50 == 0) & (ii < 270)  ) {
  currentfreq = SampledFrequency[ii]
  r = runif(1)
  if (r < currentfreq**2) {coll = 2}
  else if (r < currentfreq**2 + (1 - currentfreq)**2 ) { coll = 0}
  else {coll = 1}
  if (coll == 2 & ii > 60) {retainboolean = F; break}
  print("place1")
  
  if (coll == 0 & ii < 180)  {retainboolean = F; break}
  print("place2")
  
  if (coll == 1 & ii < 60)  {retainboolean = F; break}
  print("place3")
  
  if (coll == 1 & ii > 180) {retainboolean = F; break}
  print("place4")
  
}

}
if (retainboolean) {
  print(selcoef)
  accvals = accvals + 1
  write.table(SampledFrequency,paste0("Frequencies/Freqs", toString(accvals) ,".txt"), row.names = F, col.names = F)
}
if (accvals == TRIALS) { break}
}