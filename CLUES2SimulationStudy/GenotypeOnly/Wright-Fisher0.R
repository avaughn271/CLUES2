args = commandArgs(trailingOnly=TRUE)
N = as.numeric(args[1]) # this is haploid Ne
selectioncoef = as.numeric(args[2])
mfreq = as.numeric(args[3])
GenerationsBetweenSamples = as.numeric(args[4])

OneSim = function(modernfreq) {
  while(T) {
  TotalFreq = rep(-1, 1e7)
  TotalFreq[1] = modernfreq # This is the initial frequency
  currentfreq = TotalFreq[1]
  nextgeneration = 2
  while(currentfreq > 0 & currentfreq < 1) {
    currentfreq = rbinom(1, N, currentfreq ) / N
    TotalFreq[nextgeneration] = currentfreq
    nextgeneration = nextgeneration + 1
  }
  TotalFreq = TotalFreq[1:(nextgeneration-1)]
  if (TotalFreq[length(TotalFreq)] < 0.001) {
    return(TotalFreq)
    }
  }
}

OneTrajectory = OneSim(mfreq)
writeLines(toString(OneTrajectory[1]), "ModernFreq.txt", sep = "\n")
OneTrajectory = OneTrajectory[-1]

if (length(OneTrajectory) > 1000) {OneTrajectory = OneTrajectory[1:1000]}

M = matrix(0, nrow = 0, ncol = 4)

for (i in 1:length(OneTrajectory)) {
  if (i %% GenerationsBetweenSamples == 0) {
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
}
write.table(M, file = "Samples.txt", quote = F, sep = " ", row.names = F, col.names = F)
