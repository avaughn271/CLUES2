CHR=6

VCF = readLines(paste0("AllIndivs" ,toString(CHR), ".recode.vcf"))
NAMES = VCF[length(VCF) - 1]
SNP = VCF[length(VCF)]

NAMES = strsplit(NAMES, "\t")[[1]][-c(1,2,3,4,5,6,7,8,9)]

SNP = strsplit(SNP, "\t")[[1]][-c(1,2,3,4,5,6,7,8,9)]

PaintedNames = c()
PaintedSNP = c()
TopGenotypeLikelihood = c()
Paintings = c()
for (i in 1:length(NAMES)) {
  snpsplit = strsplit(SNP[i], ":")[[1]]
  if (snpsplit[5] != ".") {
    PaintedNames = c(PaintedNames, NAMES[i],  NAMES[i])
    PaintedSNP = c(PaintedSNP, substr(SNP[i], 1,1), substr(SNP[i], 3,3))
    Paintings = c(Paintings, substr(SNP[i],nchar(SNP[i])-2,nchar(SNP[i])-2), substr(SNP[i],nchar(SNP[i]),nchar(SNP[i])))
    
    A = strsplit(snpsplit[3], ",")[[1]]
    TopGenotypeLikelihood = c(TopGenotypeLikelihood, max(as.numeric(A[1]), as.numeric(A[2]), as.numeric(A[3])))
  }
}
hist(TopGenotypeLikelihood, main = paste0("Chromosome " ,toString(CHR)))

write.table(data.frame(PaintedNames, PaintedSNP, Paintings), quote = F, paste0("Chromosome" ,toString(CHR), ".csv"), row.names = F)