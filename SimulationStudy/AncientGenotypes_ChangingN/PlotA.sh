let generationgap=1 
let discretization=500 
let NUMTRIALS=50  ##instead of 50
mincurrfreq="0.05"
maxtimebeforepresetn=100000
let samplespersamplingtime=2

rm -f -r TrueResults
mkdir TrueResults

for ss in 0.01  0.0075 0.005  0.0025   0.001   0.0 ###for 0.0, just adjust so that simulation is back in time
do

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

if [ "0.0" = "$ss" ]; then
    maxtimebeforepresetn=1000000
fi

Rscript Wright-Fisher0.R temp $ss $generationgap $mincurrfreq $maxtimebeforepresetn $samplespersamplingtime

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --coal coal.txt --out TrueResults/output${c} --df $discretization --tCutoff 1000 --noAlleleTraj 

done

Rscript SelectionNew.R $NUMTRIALS $ss tempp $generationgap

done

Rscript LargeSelection.R $NUMTRIALS temp tempp $generationgap
mv Rplots.pdf Violin.pdf
