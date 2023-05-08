let discretization=500      ##edited this
let numimpsamples=150
totaliter=200000   #POSITION5

c=2

echo "startinfer"

time Rscript tempswap.R 1e6 7e-9 $totaliter $numimpsamples  100000  Input/samples${c} $c

echo "endinfer"

echo "startnewick"

Rscript NewicktoLengths.R $c

time python3.10 ~/desktop/CLUES2/inference.py --popFreq  0.75  --times INPUTTIMES/Times${c}.txt --N  100000  --out TrueResults/output${c} --df $discretization --tCutoff 2000 --noAlleleTraj
