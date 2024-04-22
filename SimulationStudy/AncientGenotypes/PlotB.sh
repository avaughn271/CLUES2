let N=600000 
let generationgap=1 
let discretization=2000 
let NUMTRIALS=30  ##instead of 50
mincurrfreq="0.15"
maxtimebeforepresetn=900000
let samplespersamplingtime=10

rm -f -r TrueResults
mkdir TrueResults

for ss in 0.01  0.0075 0.005  0.0025 0.001  0.00002 
do

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

if [ "0.0" = "$ss" ]; then
    maxtimebeforepresetn=10000000
fi

Rscript WrightFisherFreqMax.R $N $ss $generationgap $mincurrfreq $maxtimebeforepresetn $samplespersamplingtime

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000 --noAlleleTraj

done

Rscript SelectionNew.R $NUMTRIALS $ss $N $generationgap

done

Rscript LargeSelection.R $NUMTRIALS temp $N $generationgap
mv Rplots.pdf Violin.pdf
