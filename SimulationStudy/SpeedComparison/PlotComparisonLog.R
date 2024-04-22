OldTimes = c(7.597,  2*60 + 23.1, 
                          4*60+11.072,  
                            6*60+17.571  , 
             8*60 + 22.649 )
NewTimes = c(1.029, 5.989,10.741,13.639,17.942)
ImportanceSamples = c(1,25,50,75,100)
ymax = 600
col1 = "#e69F00"
col2 = "#009E73"

plot(ImportanceSamples, OldTimes , type="b" , bty="l" , xlab="Number of Importance Samples" , ylab="Time (seconds)" , col=col1 , lwd=3 , 
     pch=17  , ylim = c(1,1000), log='y', yaxt="n", xaxt="n")

axis(side=2, at=c(1,10,100,1000), labels = F)
text(par("usr")[1], c(1,10,100,1000),  
     labels =  c(1,10,100,1000), srt = 0, pos = 2, xpd = TRUE)

xtick<-ImportanceSamples
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xtick, srt = 0, pos = 2, xpd = TRUE)

axis(side=1, at=ImportanceSamples, labels = c("1", "25", "50", "75", "100"))

lines(ImportanceSamples , NewTimes, col=col2 , lwd=3 , pch=19 , type="b" )

# Add a legend
legend("topleft", 
       legend = c( "Without Approximations", "With Approximations"), 
       col = c(col1, col2), 
       pch = c(17,19), 
       pt.cex = 2, 
       cex = 1.2, inset = c(.01,.01))
