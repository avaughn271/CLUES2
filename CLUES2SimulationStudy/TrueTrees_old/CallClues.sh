let N=30000
let NUMTRIALS=50
modernfreq=".75"
let generationgap=1
let discretization=2000



slim BasicSelectionGenealogy.slim


python3.9 getTimesfromSLIMcap.py


python3 ~/desktop/CLUES2/inference.py --times Times.txt --popFreq 0.75 -N $N --tCutoff 1500 --df 1000 --out output


##overall, this seems to be okay. But what about the time cutoff????