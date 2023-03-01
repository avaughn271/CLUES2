ancientfreq  = runif(1, 0.05,0.2)
Generations = 500
N = 200000
selectioncoef = 0.005

M = matrix(0, nrow = 0, ncol = 4)
TotalFreq = c(ancientfreq)
currentfreq = ancientfreq
for (i in 1:Generations) {

  newmean  = (currentfreq + currentfreq * selectioncoef )  /(1 + currentfreq * selectioncoef)
  
currentfreq  = rbinom(1, N, newmean ) / N

  TotalFreq = c(TotalFreq, currentfreq)

if (i %% 1 == 0) {
    r = runif(1)
    if (r < currentfreq**2) M = rbind(M, c(paste0(toString(Generations-i), ".0001"), "-inf", "-inf", "0.00") )
    else if (r < currentfreq**2 + (1 - currentfreq)**2 ) {
      M = rbind(M, c(paste0(toString(Generations-i), ".0001"), "0.00", "-inf",  "-inf") )
    }
    else {
      M = rbind(M, c(paste0(toString(Generations-i), ".0001"), "-inf",   "0.00", "-inf") )
    }  
}}
M = M[-nrow(M),]
M = M[nrow(M):1,]
plot(TotalFreq, ylim = c(0,1))
writeLines(toString(TotalFreq ), "FrequencyTrajectory.txt", sep = "\n")
writeLines(toString(TotalFreq[length(TotalFreq)] ), "ModernFreq.txt", sep = "\n")

write.table(M, file = "Samples.txt", quote = F, sep = " ", row.names = F, col.names = F)

