library("ggplot2")

Epochs = scan("epochs.csv", sep  =",")
Midpoints = (Epochs[-1] + Epochs[-length(Epochs)])/2
Posts = read.table("logposts.csv", sep = ",")

findconfidenceinterval = function(postss) {
  actualposts = exp(postss)
  maxxindex = which(actualposts == max(actualposts))
 # print(maxxindex)
  extendd = 0
  while(TRUE) {
    #print(c(extendd, (sum(actualposts[(-extendd + maxxindex):(extendd + maxxindex)]) )  ))
    if (sum(actualposts[(-extendd + maxxindex):(extendd + maxxindex)]) > 0.95   ) {
      return(c(maxxindex-extendd, maxxindex + extendd, maxxindex))
    }
    extendd = extendd + 1
    if ((-extendd + maxxindex) < 1 | (extendd + maxxindex) > 399 ) {
      print("PROBLEM")
    }
  }
}

Lower = c()
Upper = c()
Modess = c()
for (i in 1:ncol(Posts)) {
  x = findconfidenceinterval(Posts[,i]) / nrow(Posts)
  Lower = c(Lower, x[1])
  Upper = c(Upper, x[2])
  Modess = c(Modess, x[3])
}
 


i = 1
##########################################3
  K = rev(scan(paste0("TrueResults/FrequencyTrajectory", toString(i), ".txt"), sep = ",", quiet = T))[1:length(Midpoints)]


x = rep("MAP Frequency", length(Modess))
Total = data.frame(Midpoints, Lower, Upper, Modess, x )
Modess = K
x = rep("True Frequency", length(Modess))
Total = rbind(Total, data.frame(Midpoints, Lower, Upper, Modess, x))

  ggplot(Total, aes(Midpoints, Modess, group = x)) +  geom_line(aes(color = x)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "95% Posterior Interval"), alpha = 0.3) +
    ylim(0, 1) +
    scale_shape_manual(values=c(0, .2,  1)) + 
    scale_color_manual(values=c('black','red')) +
    theme(legend.title= element_blank()) +
    scale_fill_manual("",values="black")+ labs( x="Generations Ago", y = "Derived Allele Frequency") + ggtitle("N = 30000, m = 2, s = 0.005") +
  theme(plot.title = element_text(hjust = 0.5))
  