let N=300000 
let NUMTRIALS=25
modernfreq=".75"
let generationgap=1
let discretization=2500

rm -f -r TrueResults
mkdir TrueResults

 
for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

Rscript Wright-Fisher0.R $N $modernfreq $generationgap

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000 --timeBins timebins.txt

done
 

Rscript SelectionNew.R $NUMTRIALS $ss $N $generationgap


Rscript LargeSelection.R $NUMTRIALS temp $N $generationgap
mv Rplots.pdf Violin.pdf
