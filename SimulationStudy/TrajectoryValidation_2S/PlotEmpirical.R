
TRIALS = 100
LIST = list()

for (i in 1:TRIALS) {
  LIST[[i]] = scan(paste0("Frequencies/Freqs",toString(i) , ".txt"))
}

for (frac in c(0.5,0.75,0.95,0.999))  {

Interval75Lower = rep(-1, length(LIST[[1]]))
Interval75Upper = rep(-1, length(LIST[[1]]))

for (i in 1:length(LIST[[1]])) {
  print(i)
  V = c()
  for (j in 1:length(LIST)) {
  V = c(V, LIST[[j]][i])
  }
  V = sort(V)
  currentmin = -1
  currentmax = 1000
  for (minindex in 1:(length(V) - 1)) {
    for (maxindex in (minindex+1):length(V)) {
    #  if (V[maxindex] - V[minindex] > currentmax - currentmin) {break} #added this line later
      if (sum( (V >= V[minindex])  & (V <=  V[maxindex]))  >= frac*length(V)
              & (V[maxindex] - V[minindex] < currentmax - currentmin)) {
        currentmin = V[minindex]
        currentmax = V[maxindex]
      }
    }
  }
  Interval75Lower[i] = currentmin
  Interval75Upper[i] = currentmax
}
write.table(rev(Interval75Lower), paste0("Lower", toString(frac), ".txt"), row.names = F, col.names = F  )
write.table(rev(Interval75Upper), paste0("Upper", toString(frac), ".txt")  , row.names = F, col.names = F  )
}