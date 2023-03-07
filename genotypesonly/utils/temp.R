A = scan("z_bins.txt")

B = scan("z_logcdf.txt")
C = scan("z_logsf.txt")



plot(exp(B), type = "l", col = "blue")
plot(((C)/(1-exp(B))))
