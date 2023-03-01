A = readLines("results.txt")
A = A[grep("selection", A) + 1]
library(ggplot2)
boundss = .01

coll = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",  "#CC79A7", "#D55E00", "#0072B2")


s = c()
for (i in 1:length(A)) s = c(s, as.numeric((strsplit(A[i], "\t")[[1]])[2]))

counts = rep(1,length(s) )
ss = data.frame(counts , s)
#ggplot(ss, aes(x=counts, y=s)) + geom_violin( fill = "lightskyblue" ) + 
#  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0, fill = coll[2]) + labs(  y = "Inferred #Value of s")


ggplot(ss, aes(x=counts, y=s)) + geom_violin(  fill = "lightskyblue") + 
  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0, fill = coll[2] ) +  
  ylim(-boundss, boundss) + labs(  y = "Inferred Value of s", x="") +  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()  ) + geom_hline(yintercept=0.005, linetype="dashed",  color = "red", size=1) + ggtitle("N = 200000, m = 1") +
  theme(plot.title = element_text(hjust = 0.5))

sum( (s > -boundss) * (s < boundss)) / 5000


