library(ggplot2)
par(mfrow=c(2,2))
plotlist = list()
for (sel in c(0.0025, 0.005, 0.0075, 0.01)) {
selstring = toString(sel)
CLUES2 = scan(paste0("selectioncoeff", selstring, ".txt") )
TAJIMA = scan(paste0("Tajima/Taj", selstring, ".txt") )
H12 = scan(paste0("H12/H12_", selstring, ".txt") )
nSL =scan(paste0("nSL/nSL", selstring, ".txt") )
temp = rep(0,length(nSL))
Method =  rep(c("CLUES2", "Tajima's D", "H12", "nSL"), each = length(nSL))
Values =  c(CLUES2, TAJIMA, H12, nSL)
LARGEDATAFRAME = data.frame(Values, Method)
print(sel)
plotlist[[length(plotlist) + 1]]= ggplot(LARGEDATAFRAME, aes(y=Values, x=Method , col = Method)) + 
  geom_violin(scale = "count") + 
  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0) +  
  ylim(0,0.025) +   labs(y = "Inferred Value of s", x="") + 
   ggtitle("")+ theme(legend.position="none") +
  ggtitle(paste0("s = ",selstring)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
        axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain")
        ) 
}
par(mfrow=c(2,2))

plotlist[[1]] + geom_segment(aes(x = -0.2, y = 0.0025, xend = 5, yend = 0.0025), linetype="dashed", color = "black", size=.5) 
ggsave("Final0.0025.pdf", width = 7, height = 7)

plotlist[[2]]+ geom_segment(aes(x = -0.2, y = 0.005, xend = 5, yend = 0.005), linetype="dashed", color = "black", size=.5) 
ggsave("Final0.005.pdf", width = 7, height = 7)

plotlist[[3]]+ geom_segment(aes(x = -0.2, y = 0.0075, xend = 5, yend = 0.0075), linetype="dashed", color = "black", size=.5) 
ggsave("Final0.0075.pdf", width = 7, height = 7)

plotlist[[4]]+ geom_segment(aes(x = -0.2, y = 0.01, xend = 5, yend = 0.01), linetype="dashed", color = "black", size=.5) 

ggsave("Final0.01.pdf", width = 7, height = 7)
