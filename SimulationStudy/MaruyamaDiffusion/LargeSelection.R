library(ggplot2)

ListOfSelCoeff = list()

FILENAMES = list.files(path = "./TrueResults", pattern = "selectioncoeff")
SELCOEFF = c()
for (i in 1:length(FILENAMES)) {
  ListOfSelCoeff[[i]] = scan(paste0("./TrueResults/",FILENAMES[i] ), quiet = T )
  SELCOEFF = c(SELCOEFF, substring(FILENAMES[i], 15, nchar(FILENAMES[i])  - 4   ))
}
args = commandArgs(trailingOnly=TRUE)

LARGEDATAFRAME = data.frame(matrix(nrow = 0, ncol = 2))
for (i in 1:length(ListOfSelCoeff)) {
  for (j in 1:length(ListOfSelCoeff[[1]])) {
    LARGEDATAFRAME = rbind(LARGEDATAFRAME, c(as.numeric(ListOfSelCoeff[[i]][j]),as.numeric(SELCOEFF[i] )))
  }
}
names(LARGEDATAFRAME) = c("Inferred_s", "True_s")
LARGEDATAFRAME$True_s <- as.factor(LARGEDATAFRAME$True_s)

#new
LARGEDATAFRAME$cat<-rep('True s', nrow(LARGEDATAFRAME)) 
colnames(LARGEDATAFRAME)<-c("Inferred_s","True_s","True s") #new


ggplot(LARGEDATAFRAME, aes(x=True_s, y=Inferred_s)) + geom_violin(fill = "lightskyblue", scale = "count") + 
  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0, fill = "#E69F00") + #  ylim(-0.025, 0.125) + 
 labs(y = "Inferred Value of s", x="Modern Allele Frequency") +  
  ggtitle(paste0("N = " , args[3]) )  + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_segment(aes(x = 0.6, y = 0.005, xend = 1.4, yend = 0.005, linetype=`True s`),  
                                                                color = "red", size=0.5) + 
   geom_segment(aes(x = 1.6, y = 0.005, xend = 2.4, yend =0.005, linetype=`True s`),   color = "red", size=.5) + 
  geom_segment(aes(x = 2.6, y =0.005, xend = 3.4, yend =0.005, linetype=`True s`),  color = "red", size=.5) + 
  geom_segment(aes(x = 3.6, y = 0.005, xend = 4.4, yend = 0.005, linetype=`True s`),   color = "red", size=.5) + 
  geom_segment(aes(x = 4.6, y = 0.005, xend = 5.4, yend = 0.005, linetype=`True s`),   color = "red", size=.5) +
geom_segment(aes(x = 5.6, y = 0.005, xend = 6.4, yend = 0.005, linetype=`True s`),  color = "red", size=.5) +
  theme(axis.text.x = element_text(color = "black", size = 14,   face = "plain"),
        axis.text.y = element_text(color ="black", size = 14,  face = "plain"),  
        axis.title.x = element_text(color = "black", size = 17,   face = "plain"),
        axis.title.y = element_text(color = "black", size = 17,  face = "plain"))   + #move linetype= to inside aesthetics
  scale_linetype_manual("True s",values=c("True s"=2))+
  theme(legend.title=element_blank(), 
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 4, 2, 2))