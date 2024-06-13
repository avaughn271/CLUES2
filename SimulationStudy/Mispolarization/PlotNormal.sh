 
python3.10 ~/desktop/CLUES2/inference.py --popFreq 0.65 --ancientSamps Samples.txt --N 60000 --out outputnormal --df 1500 --tCutoff 1100 --noAlleleTraj

python3.10 ~/desktop/CLUES2/inference.py --popFreq 0.35 --ancientSamps Flipped.txt --N 60000 --out outputflipped --df 1500 --tCutoff 1100 --noAlleleTraj

