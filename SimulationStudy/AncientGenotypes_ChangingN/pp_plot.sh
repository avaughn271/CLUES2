let generationgap=1 
let discretization=500
let NUMTRIALS=500  ##instead of 50
mincurrfreq="0.02"
maxtimebeforepresetn=100000
let samplespersamplingtime=2

rm -f -r TrueResults
mkdir TrueResults

ss=0

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

Rscript Wright-Fisher0.R temp 0.0 $generationgap $mincurrfreq $maxtimebeforepresetn $samplespersamplingtime

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --coal coal.txt  --out TrueResults/output${c} --df $discretization --tCutoff 1000 --noAlleleTraj

done

Rscript PP.R $NUMTRIALS $ss
mv Rplots.pdf PP.pdf

Rscript pvals.R $NUMTRIALS $ss
mv Rplots.pdf pval.pdf
