
args = commandArgs(trailingOnly=TRUE)
inputnumber = as.numeric(args[1])

options(scipen=999)
A = readLines(paste0("INPUT/Data", toString(inputnumber), ".vcf"))
if (length(A) == 6) {
  writeLines(paste0("Position", "\t", "Data"), paste0("INPUT/Data", toString(inputnumber), ".txt"))
} else {
firstline = grep("CHROM", A) + 1
numberhaplotypes = (length(strsplit(A[firstline], "\t")[[1]]) - 9) * 2
NAMESLINE = paste0("Position", "\t", "Data")

SITEVECTOR = rep("temp", length(A) - firstline + 1)

for (i in firstline:length(A)) {
  totalstring = ""
  SPLIT = strsplit(A[i], "\t")[[1]]
  position = SPLIT[2]
  for (j in 10:length(SPLIT)) {
    if (SPLIT[j] == "0|0") {totalstring = paste0(totalstring, "00")}
    else if (SPLIT[j] == "1|0") {totalstring = paste0(totalstring , "10")}
    else if (SPLIT[j] == "0|1") {totalstring = paste0(totalstring , "01")}
    else {totalstring = paste0(totalstring ,"11")}
  }
  SITEVECTOR[i - firstline + 1] = paste0(position, "\t", totalstring)
}

writeLines(c(NAMESLINE, SITEVECTOR), paste0("INPUT/Data", toString(inputnumber), ".txt"))

}