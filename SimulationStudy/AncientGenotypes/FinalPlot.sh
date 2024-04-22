let N=10000 
let NUMTRIALS=200
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

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 500

done

Rscript PP.R $NUMTRIALS $ss
mv Rplots.pdf PP.pdf

Rscript pvals.R $NUMTRIALS $ss
mv Rplots.pdf pval.pdf


##compare with null with original code, which should work. Run it on the exact same data and see what the difference is, if there is one. Options are:

#improper optimization algorithm
#tcutoff is different
#data generation process itself, possible if the answers they give are quite similar.
#probability cutoffs.
# exit probabilities are different.
#