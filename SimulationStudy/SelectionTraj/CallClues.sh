
rm -f results.txt
touch results.txt

rm -f -r TrueResults
mkdir TrueResults

for (( c=1; c<=1; c++ ))
do

echo $c

Rscript Wright-Fisher.R

mv Samples.txt TrueResults/Samples${c}.txt
mv FrequencyTrajectory.txt TrueResults/FrequencyTrajectory${c}.txt
mv ModernFreq.txt TrueResults/ModernFreq${c}.txt

line=$(head -n 1 TrueResults/ModernFreq${c}.txt)

python3 clues-master/inference.py --popFreq ${line}  --ancientSamps TrueResults/Samples${c}.txt --N 30000 --df 400  --out output

done

mv Rplots.pdf truth.pdf

python3 plot_traj.py


Rscript PlotEnvelope.R
