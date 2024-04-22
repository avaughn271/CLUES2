let N=10000 
let NUMTRIALS=450
let generationgap=2
let discretization=450
mincurrfreq="0.02"
maxtimebeforepresetn=100000
let samplespersamplingtime=1

rm -f -r TrueResults
mkdir TrueResults

ss=0

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

Rscript Wright-Fisher.R $N 0.0 $generationgap $mincurrfreq $maxtimebeforepresetn $samplespersamplingtime

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000 --timeBins 200 600 --noAlleleTraj

done

Rscript PP.R $NUMTRIALS $ss
mv Rplots.pdf PP.pdf

Rscript pvals.R $NUMTRIALS $ss
mv Rplots.pdf pval.pdf