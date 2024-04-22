CHR = 16

Data = read.table(paste0("Chromosome", toString(CHR), ".csv"), header = T)
INFO = read.csv("../StoneAgeData/MesoNeoData.tsv", sep = "\t")
Dates = c()

for (i in 1:nrow(Data)) {
  samplename = Data$PaintedNames[i]
  Dates = c(Dates, INFO$ageAverage[which(INFO$sampleId == samplename)] / 28.0)
}

Data$Age = Dates

TimeOrder = order(Dates)

Data = Data[TimeOrder, ]

if (CHR == 16) {Data$PaintedSNP = 1 - Data$PaintedSNP}
ToWrite = ""

for (i in 1:nrow(Data)) {
  if (Data$PaintedSNP[i] == 1) {
    ToWrite = paste0(ToWrite, toString(Data$Age[i]),  " -inf 0", "\n")
  } else if (Data$PaintedSNP[i] == 0) {
    ToWrite = paste0(ToWrite,  toString(Data$Age[i]),  " 0 -inf", "\n")
  }
  else {print("PROBLEM")}
}
writeLines(ToWrite, paste0("Haplotypes", toString(CHR), ".csv"))
