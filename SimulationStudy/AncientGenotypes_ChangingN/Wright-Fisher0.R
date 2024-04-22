recentpopsize = 60000
middlepopsize = 200000
oldpopsize = 80000
#400 and 600 are the changepoints

args = commandArgs(trailingOnly=TRUE)
temp = (args[1]) # this is haploid Ne
selectioncoef = as.numeric(args[2])
h = 0.5 ##1, 1 + hs, 1 + s
GenerationsBetweenSamples = as.numeric(args[3])
mincurrfrequnecy = as.numeric(args[4])

MAXTIME_BeforePresent =  as.numeric(args[5])
samplespersamplingtime =  as.numeric(args[6])

OneTrajectoryFunction = function() {
  while(T) {
  timeofmutation = sample(1:MAXTIME_BeforePresent, 1)
  if (timeofmutation < 400 & runif(1) > 3/10) {next} # mutation happens proportional to popsize
  if (timeofmutation > 600 & runif(1) > 4/10) {next}
  TotalFreq = rep(-1, timeofmutation)
  if (timeofmutation < 400) {TotalFreq[1] = 1/recentpopsize }
  if (timeofmutation >= 400 & timeofmutation < 600) {TotalFreq[1] = 1/middlepopsize}
  if (timeofmutation >= 600) {TotalFreq[1] = 1/oldpopsize}
  currentfreq = TotalFreq[1]
  nextgeneration = 2
  while(currentfreq > 0 & currentfreq < 1) {
if (nextgeneration > (timeofmutation -  400    )  )  {N = 60000} else 
if (nextgeneration > (timeofmutation -  600    )      )  {N = 200000} else {
     N = 80000   }
#print(c(N, nextgeneration - timeofmutation ))
#print(c(N, nextgeneration, timeofmutation))
    p = currentfreq
    newmean = p - (selectioncoef* (-1 + p)* p* (-p + h *(-1 + 2* p))) /(-1 + selectioncoef *(2* h *(-1 + p) - p) *p)
    currentfreq = rbinom(1, N, newmean) / N
    TotalFreq[nextgeneration] = currentfreq
    nextgeneration = nextgeneration + 1
    if (   nextgeneration > (length(TotalFreq) + 1) ) {next}
    if (  (nextgeneration == (length(TotalFreq) + 1) )  &  (currentfreq >= mincurrfrequnecy) ){ # enforce that modern freq is at least 25%.
#print(c(timeofmutation, currentfreq, length(TotalFreq)))
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