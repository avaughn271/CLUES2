let N=40000      #POSITION1
let NUMTRIALS=30  #POSITION2
modernfreq=".75"
let discretization=450      ##edited this
let numleaves=24    #changed this tooooo
LL="1e6"
muu="4e-9"       #number1
let numimpsamples=600
 

for ss in   0.0 0.001 0.01 0.0025 0.005 0.0075
do


rm -f -r INPUT
mkdir INPUT

python3.9 GetData.py $N $numleaves $ss $NUMTRIALS $LL $muu
mv INPUT INPUT${ss}
done
