let N=40000
let NUMTRIALS=50 #30
let discretization=450
let numleaves=20
let tcuto=3000

rm -f -r TrueResults
mkdir TrueResults

for ss in 0.01 0.001 0.0025 0.005 0.0075  0.0
do

rm -f -r INPUTTIMES
mkdir INPUTTIMES


for (( c=1; c<=${NUMTRIALS}; c++ ))
do

Rscript Emissions.R $N $ss temp 0.1  300000 $numleaves

echo $c
line=$(head -n 1 ModernFreq.txt)

mv INPUTTIMES/AncientTimes.txt  INPUTTIMES/AncientTimes${c}.txt
mv INPUTTIMES/Samples.txt  INPUTTIMES/Samples${c}.txt

python3.10 ~/desktop/CLUES2/inference.py --popFreq $line --times INPUTTIMES/AncientTimes${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff $tcuto --noAlleleTraj --ancientSamps INPUTTIMES/Samples${c}.txt

done

Rscript SelectionNew.R $NUMTRIALS $ss

done

Rscript LargeSelection.R $NUMTRIALS temp $N
mv Rplots.pdf Violin.pdf