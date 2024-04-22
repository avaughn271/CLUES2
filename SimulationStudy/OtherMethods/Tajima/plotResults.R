A = scan("selcoeffs.txt")

denss  = density(A, kernel="rectangular", bw = 0.005)
plot(denss)
denss$x[which(denss$y == max(denss$y))]


denss  = density(A,  bw = 0.002)
plot(denss)
denss$x[which(denss$y == max(denss$y))]


denss  = density(A)
plot(denss)
denss$x[which(denss$y == max(denss$y))]

hist(A)