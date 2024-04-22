let N=40000
let NUMTRIALS=30  #POSITION1
modernfreq=".75"
let numleaves=24
LL="1e6"
muu="4e-9"
numsamples=1000 ##100,0.1 is 3 hours
epsilon=0.05

rm -f -r TrueResults
mkdir TrueResults
 
for ss in 0.005 0.0025# 0.01 0.0075
do

python3.9 getnSL.py $ss $numleaves

max_parallel_processes=30



for (( c=0; c<${NUMTRIALS}; c++ ))
do 
{
taj=$(head -n 1 TrueResults/nSL_${c}.txt)

echo $taj
python3.9 RejectionnSL.py $N $numleaves $ss $NUMTRIALS $LL $muu $numsamples $taj $epsilon $c
} &

done

wait
Rscript plotAllResultsnSL.R $ss

done


Rscript DEBUGnSL.R