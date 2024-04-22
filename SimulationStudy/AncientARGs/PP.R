args = commandArgs(trailingOnly=TRUE)
NumberOfTrials = as.numeric(args[1])
TrueS = args[2]
library("ggplot2")

loglr = c()
for (i in 1:NumberOfTrials) {
  Line1 = readLines(paste0("TrueResults/output",toString(i),"_inference.txt"))[1]
  
  loglr = c(loglr, as.numeric(read.table(paste0("TrueResults/output",toString(i),"_inference.txt" ))[[1]][2])   )
}
  

vals = 2*loglr
desiredquantiles = c(qchisq(c(1:1000)/1000, 1) , 0, min(vals)-0.1, max(vals)+0.1, ( (sort(vals)[-1] + sort(vals)[-length(vals)])/2.0) )

actual = c()
experimental = c()
for (i in desiredquantiles) {
  actual = c(actual, pchisq(i, 1))
  experimental = c(experimental, mean(vals <= i))
}

plot(actual,experimental, xlim = c(0,1), ylim = c(0,1), pch = 19, col = "steelblue4",
     main = "", xlab = "Theoretical Distribution Function", ylab = "Empirical Distribution Function")
lines(c(0,1), c(0,1), col = "red", lty=2,lwd=2)
