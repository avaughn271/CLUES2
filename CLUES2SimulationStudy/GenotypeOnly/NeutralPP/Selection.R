library(ggplot2)
boundss = .01

args = commandArgs(trailingOnly=TRUE)

s = c()
for (i in 1:as.numeric(args[1])) {
  Line3 = readLines(paste0("TrueResults/output",toString(i),"_inference.txt"))[3]
  s = c(s, as.numeric((strsplit(Line3, "\t")[[1]])[2]))
}

counts = rep(1,length(s))
ss = data.frame(counts, s)

ggplot(ss, aes(x=counts, y=s)) + geom_violin(fill = "lightskyblue") + 
  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0, fill = "#E69F00") +  
  ylim(min(c(-boundss, s-0.0001)), max(c(boundss,s+0.0001))) + labs(y = "Inferred Value of s", x="") +  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())  +
  geom_hline(yintercept=as.numeric(args[2]), linetype="dashed", color = "red", size=1) +
  ggtitle(paste0("N = ",toString(args[3]), ", m = ",toString(args[4]))) + 
  theme(plot.title = element_text(hjust = 0.5))
