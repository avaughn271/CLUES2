let N=30000
let NUMTRIALS=10
modernfreq=".75"
let discretization=600
let numleaves=8

rm -f -r TrueResults
mkdir TrueResults

var2=2
DipNe=$((N / var2))


for ss in 0.0 0.01 0.001 0.0025 0.005 0.0075  ##0.0 have to do separately.
do

rm -f -r INPUTTIMES
mkdir INPUTTIMES

python3.9 recapitateall.py $N $numleaves $ss $NUMTRIALS

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

python3.10 ~/desktop/CLUES2/inference.py --popFreq  $modernfreq --times INPUTTIMES/Times${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000

done

Rscript Selection.R $NUMTRIALS $ss

done

Rscript LargeSelection.R $NUMTRIALS temp $N
mv Rplots.pdf Violin.pdf