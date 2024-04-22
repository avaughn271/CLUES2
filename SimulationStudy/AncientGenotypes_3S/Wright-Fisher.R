args = commandArgs(trailingOnly=TRUE)
N = as.numeric(args[1]) # this is haploid Ne
selectioncoef = as.numeric(args[2])
h = 0.5 ##1, 1 + hs, 1 + s
GenerationsBetweenSamples = as.numeric(args[3])
mincurrfrequnecy = as.numeric(args[4])

MAXTIME_BeforePresent =  as.numeric(args[5])
samplespersamplingtime =  as.numeric(args[6])

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
    if (  (nextgeneration > length(TotalFreq)) &    (currentfreq > mincurrfrequnecy) ){ # enforce that modern freq is at least 25%.
print(c(timeofmutation, currentfreq))
      return(TotalFreq)
    }
  }
  }
}

OneTrajectory = rev(OneTrajectoryFunction())
#print(OneTrajectory)
writeLines(toString(OneTrajectory[1]), "ModernFreq.txt", sep = "\n")
OneTrajectory = OneTrajectory[-1]
if (length(OneTrajectory) < 1000) {OneTrajectory = c(OneTrajectory, rep(0,10000))}
if (length(OneTrajectory) > 1000) {OneTrajectory = OneTrajectory[1:1000]}
M = matrix(0, nrow = 0, ncol = 4)

for (i in 1:length(OneTrajectory)) {
  if (i %% GenerationsBetweenSamples == 0) {
    for (tempindex in 1:samplespersamplingtime) {
    freq = OneTrajectory[i]
      r = runif(1)
      if (r < freq**2) M = rbind(M, c(paste0(toString(i), ".0001"), "-inf", "-inf", "0.00") )
      else if (r < freq**2 + (1 - freq)**2 ) {
        M = rbind(M, c(paste0(toString(i), ".0001"), "0.00", "-inf",  "-inf") )
      }
      else {
        M = rbind(M, c(paste0(toString(i), ".0001"), "-inf",   "0.00", "-inf") )
      }
      }
}}
write.table(M, file = "Samples.txt", quote = F, sep = " ", row.names = F, col.names = F)