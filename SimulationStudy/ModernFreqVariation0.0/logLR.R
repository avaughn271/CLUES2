args = commandArgs(trailingOnly=TRUE)
s = c()
for (i in 1:as.numeric(args[1])) {
  Line3 = readLines(paste0("TrueResults/output",toString(i),"_inference.txt"))[2]
  s = c(s,((strsplit(Line3, "\t")[[1]])[1]))
}

writeLines(s, paste0("TrueResults/selectioncoeff", toString(args[2])  ,".txt"))
