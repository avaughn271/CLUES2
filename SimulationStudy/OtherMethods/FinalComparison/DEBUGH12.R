library(ggplot2)
TAJIMA1 = scan("H12_0.0025.txt")
TAJIMA2 = scan("H12_0.005.txt")
TAJIMA3 = scan("H12_0.0075.txt")
TAJIMA4 = scan("H12_0.01.txt")

temp = rep(0,length(TAJIMA1))
Method =  rep(c("Taj1", "Taj2", "Taj3", "Taj4"), each = length(TAJIMA1))
Values =  c(TAJIMA1,TAJIMA2,TAJIMA3,TAJIMA4)
LARGEDATAFRAME = data.frame(Values, Method)

ggplot(LARGEDATAFRAME, aes(y=Values, x=Method , col = Method)) + geom_violin(scale = "count") + 
  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0) +  
  ylim(0,0.025) +   labs(y = "Inferred Value of s", x="") + 
  ggtitle("")+ theme(legend.position="none")

##100,0.5, 1 minute
