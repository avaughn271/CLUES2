let NUMTRIALS=50
ss=0.005
let discretization=900
let numleaves=200

for N in 100000 10000 1000 100
do

rm -f -r TrueResults
mkdir TrueResults

for modernfreq in 0.99  0.95 0.7 0.3 0.05 0.01 
do

rm -f -r INPUTTIMES
mkdir INPUTTIMES

python3.9 recapitateallnew.py $N $numleaves $ss $NUMTRIALS $modernfreq

for (( c=1; c<=${NUMTRIALS}; c++ ))
do

python3.10 ~/desktop/CLUES2/inference.py --popFreq  $modernfreq --times INPUTTIMES/Times${c}.txt --N $N --out TrueResults/output${c} --df $discretization --tCutoff 10000 --noAlleleTraj #changed cutoff from 1000 to 10000

done

Rscript logLR.R $NUMTRIALS $modernfreq

done

#Rscript LargeSelection.R $NUMTRIALS temp $N
#mv Rplots.pdf Violin${N}.pdf
mv TrueResults TrueResultslogLR${N}

done




for N in 100000 10000 1000 100
do
 
Rscript plotlogLR.R $N
 mv Rplots.pdf logLR${N}.pdf
done

