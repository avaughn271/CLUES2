let N=30000 
let NUMTRIALS=50
modernfreq=".75"
let generationgap=4  
let discretization=450 

rm -f -r TrueResults
mkdir TrueResults

for ss in 0.0 0.01 0.001 0.0025 0.005 0.0075 ###for 0.0, just adjust so that simulation is back in time
do


for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c



if [ "0.0" = "$ss" ]; then
    Rscript Wright-Fisher0.R $N $ss $modernfreq $generationgap
else
      Rscript Wright-Fisher.R $N $ss $modernfreq $generationgap
fi


mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3.10 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 1000

done

#Rscript PP.R $NUMTRIALS $ss
#mv Rplots.pdf PP${ss}.pdf

Rscript SelectionNew.R $NUMTRIALS $ss $N $generationgap

done

Rscript LargeSelection.R $NUMTRIALS temp $N $generationgap
mv Rplots.pdf Violin.pdf
