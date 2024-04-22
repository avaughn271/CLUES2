Approx = read.table("approx.txt")
Exact = read.table("exact.txt")
width = 5
col2 = "#e69F00"
col1 = "#009E73"

plot(Approx[[1]], Approx[[2]],  xlab="s" , ylab="log-likelihood", 
     col=col1 , lwd=width, type = "l"  , ylim = c(1,7), xaxt = "n", xlim = c(0.01,0.041))
axis(1, at = c(seq(.012,0.04, by = 0.004)))
lines(Exact[[1]] , Exact[[2]], col=col2 , lwd=width,lty = 2)

# Add a legend
legend("topleft", 
       legend =rev( c("With Approximations", "Without Approximations")), 
       col = rev( c(col1, col2)), 
       lty = rev( c(1,2)), lwd = 3, 
       cex = 1, inset = c(.01,.01))
