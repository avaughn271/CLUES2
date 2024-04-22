let N=50000  # moved from 50000 to 10000
let generationgap=1 
let discretization=500 
let NUMTRIALS=50  ##instead of 50
mincurrfreq="0.3" #and set max to be 0.7
maxtimebeforepresetn=100000 #changed frfom 100000
let samplespersamplingtime=8

rm -f -r TrueResults
mkdir TrueResults

for ss in     0.0 1.0 2.0
do

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

if [ "2.0" = "$ss" ]; then
    maxtimebeforepresetn=1000
fi



Rscript WrightFisherFreqMax.R $N $ss $generationgap $mincurrfreq $maxtimebeforepresetn $samplespersamplingtime

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000 --noAlleleTraj --h $ss

done

Rscript SelectionNew.R $NUMTRIALS $ss $N $generationgap

done

Rscript LargeSelection.R $NUMTRIALS temp $N $generationgap
mv Rplots.pdf Violin.pdf
