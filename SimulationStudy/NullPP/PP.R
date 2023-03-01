A = readLines("results.txt")
A = A[grep("log", A)]
library("ggplot2")

loglr = c()
for (i in 1:length(A)) loglr = c(loglr, as.numeric(substr(A[i], 8, nchar(A[i]) )))

vals = 2*loglr


desiredquantiles = c(qchisq(c(1:10000)/10000, 1) , 0,  seq(min(vals), max(vals), length.out = 10000))


actual = c()
experimental = c()
for (i in desiredquantiles) {
  actual = c(actual, pchisq(i , 1))
  experimental = c(experimental, mean(vals <= i))
}

plot(actual,experimental, xlim = c(0,1), ylim = c(0,1), pch = 19, col = "steelblue4", main = "s = 0", xlab = "Theoretical Distribution Function", ylab = "Empirical Distribution Function")
lines(c(-100,100), c(-100,100), col = "red", lty=2)
