args = commandArgs(trailingOnly=TRUE)
NumberOfTrials = as.numeric(args[1])
TrueS = args[2]
library("ggplot2")

pvals = c()
for (i in 1:NumberOfTrials) {
  Line1 = readLines(paste0("TrueResults/output",toString(i),"_inference.txt"))[1]
  
  pvals = c(pvals, 10**(-as.numeric(read.table(paste0("TrueResults/output",toString(i),"_inference.txt" ))[[2]][2])   ))
}

hist(pvals, breaks = 5, xlab = "p-values", main = "", col = "lightskyblue")
lines(c(-10000,1000), c(NumberOfTrials/5, NumberOfTrials/5) , lty =2, col = "red")