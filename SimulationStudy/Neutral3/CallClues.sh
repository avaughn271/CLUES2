
rm -f results.txt
touch results.txt

rm -f -r TrueResults
mkdir TrueResults

for (( c=1; c<=50; c++ ))
do

echo $c

Rscript Wright-Fisher.R

mv Samples.txt TrueResults/Samples${c}.txt
mv FrequencyTrajectory.txt TrueResults/FrequencyTrajectory${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

python3 clues-master/inference.py --popFreq ${line}  --ancientSamps TrueResults/Samples${c}.txt --N 200000 --df 600 >> results.txt

done

Rscript Selection.R
mv Rplots.pdf selection.pdf

