let N=100000      #POSITION1
let NUMTRIALS=30  #POSITION2
modernfreq=".75"
let discretization=1500      ##edited this
let numleaves=40    #changed this tooooo
LL="1e6"              
muu="7e-7"       #number1
let numimpsamples=1
totaliter=30000   #changed!!!

rm -f -r TrueResults
mkdir TrueResults

var2=2
DipNe=$((N / var2))


rm -f -r INPUT
mkdir INPUT

for ss in  0.01 0.0 0.001 0.0025 0.005 0.0075
do

rm -f -r INPUTTIMES
mkdir INPUTTIMES

python3.9 GetData.py $N $numleaves $ss $NUMTRIALS $LL $muu

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

time python3.10 ~/desktop/CLUES2/inference.py --popFreq  $modernfreq --times INPUTTIMES/Times${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 6000   #CHANGEEEE1
echo "endclues"

done

Rscript Selection.R $NUMTRIALS $ss

done

Rscript LargeSelection.R $NUMTRIALS temp $N
mv Rplots.pdf Violin.pdf
