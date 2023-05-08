args = commandArgs(trailingOnly=TRUE)
recentsel = c()
middlesel = c()
oldsel = c()
for (i in 1:as.numeric(args[1])) {
  Line3 = readLines(paste0("TrueResults/output",toString(i),"_inference.txt"))[2]
  recentsel = c(recentsel, strsplit(Line3, "\t")[[1]][5])
middlesel = c(middlesel,  strsplit(Line3, "\t")[[1]][8])
oldsel = c(oldsel,  strsplit(Line3, "\t")[[1]][11])

}

writeLines(recentsel, paste0("TrueResults/selectioncoeff0.01.txt"))
writeLines(middlesel, paste0("TrueResults/selectioncoeff0.0.txt"))
writeLines(oldsel, paste0("TrueResults/selectioncoeff-0.005.txt"))
