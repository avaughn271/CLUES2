let N=40000 
let NUMTRIALS=600 #increase later
let discretization=450
mincurrfreq="0.05"
let numleaves=20
let tcuto=1000

rm -f -r TrueResults
mkdir TrueResults

ss=0


rm -f -r INPUTTIMES
mkdir INPUTTIMES

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

Rscript AncientARGs.R $N $ss temp  $mincurrfreq  300000 $numleaves

echo $c
line=$(head -n 1 ModernFreq.txt)
mv INPUTTIMES/AncientTimes.txt  INPUTTIMES/AncientTimes${c}.txt


python3.10 ~/desktop/CLUES2/inference.py --popFreq  $line --times INPUTTIMES/AncientTimes${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff  $tcuto --noAlleleTraj

done

Rscript PP.R $NUMTRIALS $ss
mv Rplots.pdf PP.pdf

Rscript pvals.R $NUMTRIALS $ss
mv Rplots.pdf pval.pdf
