let N=30000 
let NUMTRIALS=25
modernfreq=".75"
let generationgap=1
let discretization=2000 

rm -f -r TrueResults
mkdir TrueResults

####ESTIMATES_ARE_GODD_CHECK_LogLR, if not, replace coal time with something else and find problem!!

for ss in 0.0 0.01 0.001 0.0025 0.005 0.0075
do


for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c


Rscript Wright-Fisher0.R $N $ss $modernfreq $generationgap

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt  --coal coal.txt --out TrueResults/output${c} --df $discretization --tCutoff 1000

done

Rscript SelectionNew.R $NUMTRIALS $ss $N $generationgap

done

Rscript LargeSelection.R $NUMTRIALS temp $N $generationgap
mv Rplots.pdf Violin.pdf
