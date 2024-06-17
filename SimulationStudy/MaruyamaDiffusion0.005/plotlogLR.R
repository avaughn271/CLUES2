library(ggplot2)
#3 locations to change for each population
ListOfSelCoeff = list()
args = commandArgs(trailingOnly=TRUE)
foldername = paste0("./TrueResultslogLR", args[1])
FILENAMES = list.files(path = foldername, pattern = "selectioncoeff")
SELCOEFF = c()
for (i in 1:length(FILENAMES)) {
  ListOfSelCoeff[[i]] = scan(paste0(foldername, "/", FILENAMES[i] ), quiet = T )
  SELCOEFF = c(SELCOEFF, substring(FILENAMES[i], 15, nchar(FILENAMES[i])  - 4   ))
}

LARGEDATAFRAME = data.frame(matrix(nrow = 0, ncol = 2))
for (i in 1:length(ListOfSelCoeff)) {
  for (j in 1:length(ListOfSelCoeff[[1]])) {
    temp  = as.numeric(ListOfSelCoeff[[i]][j])
   # LARGEDATAFRAME = rbind(LARGEDATAFRAME, c(  1-pchisq(  temp + temp, 1 )    ,as.numeric(SELCOEFF[i] )))
    LARGEDATAFRAME = rbind(LARGEDATAFRAME, c(  temp  , as.numeric(SELCOEFF[i] )))
    
     }
}
names(LARGEDATAFRAME) = c("Inferred_s", "True_s")
LARGEDATAFRAME$True_s <- as.factor(LARGEDATAFRAME$True_s)

#new
LARGEDATAFRAME$cat<-rep('True s', nrow(LARGEDATAFRAME)) 
colnames(LARGEDATAFRAME)<-c("Inferred_s","True_s","True s") #new


ggplot(LARGEDATAFRAME, aes(x=True_s, y=Inferred_s)) + geom_violin(fill = "lightskyblue", scale = "count") + 
  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0, fill = "#E69F00") +   ylim(-2, 10.0) + 
 labs(y = "Log-likelihood Ratio", x="Modern Allele Frequency") +  
  ggtitle(paste0("N = " ,  args[1] ) )  +  
  theme(plot.title = element_text(hjust = 0.5)) +  
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