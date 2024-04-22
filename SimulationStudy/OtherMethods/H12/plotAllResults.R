B = c()
for (i in 1:30) {
A = scan(paste0("TrueResults/selcoeffs", toString(i-1) , ".txt" ))

denss = density(A)

B = c(B, denss$x[which(denss$y == max(denss$y))] )

}

trueselection = as.numeric(commandArgs(trailingOnly = TRUE))[1]

boxplot(B)
print(B)
lines(c(-1000,10000), c(trueselection,  trueselection))