let N=30000
let NUMTRIALS=50
ss=".005"      ###should simulate   0, 0.001, 0.0025, 0.005, 0.0075, 0.01
modernfreq=".75"
let generationgap=4


rm -f -r TrueResults
mkdir TrueResults

for ss in 0.01 0.001 0.0025 0.005 0.0075 ###for 0.0, just adjust so that simulation is back in time
do


for (( c=1; c<=${NUMTRIALS}; c++ ))
do

echo $c

Rscript Wright-Fisher.R $N $ss $modernfreq $generationgap

mv Samples.txt TrueResults/Samples${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

echo "runningclues"

python3 ~/desktop/CLUES2/inference.py --popFreq ${line} --ancientSamps TrueResults/Samples${c}.txt --N $N --out TrueResults/output${c} --df 450

done

Rscript PP.R $NUMTRIALS $ss
mv Rplots.pdf PP${ss}.pdf

Rscript Selection.R $NUMTRIALS $ss $N $generationgap

mv Rplots.pdf Violin${ss}.pdf


done