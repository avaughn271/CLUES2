library(ggplot2)
boundss = .01

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
ggplot(LARGEDATAFRAME, aes(x=True_s, y=Inferred_s)) + geom_violin(fill = "lightskyblue", scale = "count") + 
  geom_boxplot(width=0.05, outlier.shape = NA, coef = 0, fill = "#E69F00") +  
  ylim(min(c(-boundss, LARGEDATAFRAME[[1]]-0.0001)), max(c(boundss,LARGEDATAFRAME[[1]]+0.0001))) + labs(y = "Inferred Value of s", x="True Values of s") +  
  ggtitle(paste0("N = ",toString(args[3]), ", m = ",toString(args[4]))) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_segment(aes(x = 0.6, y = -0.005, xend = 1.4, yend =-0.005), linetype="dashed", color = "red", size=0.5) +
   geom_segment(aes(x = 1.6, y =  0.0, xend = 2.4, yend = 0.0), linetype="dashed", color = "red", size=.5) + 
  geom_segment(aes(x = 2.6, y = 0.01, xend = 3.4, yend = 0.01), linetype="dashed", color = "red", size=.5) 
