args = commandArgs(trailingOnly=TRUE)

s = c()
for (i in 1:as.numeric(args[1])) {
  Line3 = readLines(paste0("TrueResults/output",toString(i),"_inference.txt"))[3]
  s = c(s,((strsplit(Line3, "\t")[[1]])[2]))
}

writeLines(s, paste0("TrueResults/selectioncoeff", toString(args[2])  ,".txt"))
