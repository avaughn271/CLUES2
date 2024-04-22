let N=30000 
let NUMTRIALS=1000
modernfreq=".75"
let numleaves=80
let discretization=650

rm -f -r TrueResults
mkdir TrueResults

ss=0

rm -f -r INPUTTIMES
mkdir INPUTTIMES

python3.9 recapitateallnew.py $N $numleaves $ss $NUMTRIALS

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq $modernfreq --times  INPUTTIMES/Times${c}.txt   --N $N --out TrueResults/output${c} --df $discretization --tCutoff 500 --noAlleleTraj

done

Rscript PP.R $NUMTRIALS $ss
mv Rplots.pdf PP.pdf

Rscript pvals.R $NUMTRIALS $ss
mv Rplots.pdf pval.pdf
