let N=70000 
let NUMTRIALS=30
modernfreq=".75"
let generationgap=1
let discretization=600 ##2500
mincurrfreq="0.1"
maxtimebeforepresetn=100000
let samplespersamplingtime=8

rm -f -r TrueResults
mkdir TrueResults

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

Rscript Wright-Fisher0.R $N 0.01 $generationgap $mincurrfreq $maxtimebeforepresetn  $samplespersamplingtime


mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000 --timeBins  200 600 --noAlleleTraj

done
 
Rscript SelectionNew.R $NUMTRIALS $ss $N $generationgap

Rscript LargeSelection.R
mv Rplots.pdf Violin.pdf
