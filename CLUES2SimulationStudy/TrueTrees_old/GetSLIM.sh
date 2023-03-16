let N=30000
let NUMTRIALS=50
modernfreq=".75"
let generationgap=1
let discretization=600




rm -r -f SLIMTREES
mkdir SLIMTREES


x=1
while [ $x -lt $NUMTRIALS ]
do


slim treeseq.slim

  x=($(ls SLIMTREES | wc -l))
echo $x
done


rm -f -r INPUTTIMES
mkdir INPUTTIMES


rm -r  -f TrueResults
mkdir TrueResults

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

python3.10 ~/desktop/CLUES2/inference.py --popFreq  $modernfreq --times INPUTTIMES/Times${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000

done



Rscript Selection.R $NUMTRIALS 0.01


Rscript LargeSelection.R $NUMTRIALS temp $N
mv Rplots.pdf Violin.pdf