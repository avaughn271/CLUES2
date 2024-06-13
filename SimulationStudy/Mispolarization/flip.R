A = read.table("Samples.txt")

head(A)

B = cbind(A[[1]], A[[4]], A[[3]], A[[2]])

write.table(B, "Flipped.txt", quote = F,
            row.names= F, col.names = F)