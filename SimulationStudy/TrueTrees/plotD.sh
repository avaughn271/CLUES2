let N=600000
let NUMTRIALS=25
modernfreq=".75"
let discretization=2500
let numleaves=800

rm -f -r TrueResults
mkdir TrueResults

var2=2
DipNe=$((N / var2))


for ss in 0.01 0.001 0.0025 0.005 0.0075  0.0 ##0.0 have to do separately.
do

rm -f -r INPUTTIMES
mkdir INPUTTIMES

python3.9 recapitateallnew.py $N $numleaves $ss $NUMTRIALS

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

python3.10 ~/desktop/CLUES2/inference.py --popFreq  $modernfreq --times INPUTTIMES/Times${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000

done

Rscript SelectionNew.R $NUMTRIALS $ss

done

Rscript LargeSelection.R $NUMTRIALS temp $N
mv Rplots.pdf Violin.pdf