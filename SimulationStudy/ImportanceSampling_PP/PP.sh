let N=20000
let NUMTRIALS=1000  #POSITION1
modernfreq=".5"
let discretization=450
let numleaves=12
LL="1e6"
muu="7e-9"
let numimpsamples=100 #POSITION2 # cahnged from 75
totaliter=75000  #POSITION3 # changed from 7500

rm -f -r TrueResults
mkdir TrueResults

rm -f -r INPUT
mkdir INPUT

rm -f -r INPUTTIMES
mkdir INPUTTIMES

python3.9 GetData.py $N $numleaves 0.0 $NUMTRIALS $LL $muu

for (( c=1; c<=$NUMTRIALS; c++ ))
do

Rscript VCFtoInput.R $c

rm -f -r Input/samples${c}
mkdir  Input/samples${c}

echo "startinfer"

time Rscript tempswap.R $LL $muu $totaliter $numimpsamples $N Input/samples${c} $c

echo "endinfer"

echo "startnewick"

Rscript NewicktoLengths.R $c

echo "startclues"

time python3.10 ~/desktop/CLUES2/inference.py --popFreq  $modernfreq --times INPUTTIMES/Times${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000  --noAlleleTraj
echo "endclues"

done

Rscript PP.R $NUMTRIALS $ss
mv Rplots.pdf PP.pdf

Rscript pvals.R $NUMTRIALS $ss
mv Rplots.pdf pval.pdf
