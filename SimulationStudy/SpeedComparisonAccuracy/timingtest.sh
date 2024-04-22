
time python3.10 inference_approx_new.py       --popFreq 0.75 --times Times100.txt --N 80000 --out CLUESapprox --df 450 --tCutoff 400 > approx.txt

time python3.10 inference_exact.py            --popFreq 0.75 --times Times100.txt --N 80000 --out CLUESapprox --df 450 --tCutoff 400 > exact.txt

Rscript PlotComparison.R