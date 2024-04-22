
python3.10 ~/desktop/CLUES2/inference.py --popFreq 0.75 --ancientSamps TrueSamplesNew.txt --N 30000 --out output --tCutoff 403 --df 600 --integration_points 100

python3.10 ~/desktop/CLUES2/plot_traj.py --freqs output_freqs.txt --post output_post.txt --figure output --posterior_intervals 0.5 0.75 0.95

rm -f -r Frequencies

mkdir Frequencies

time Rscript Emissions_old.R

Rscript PlotEmpiricalNew.R

python3.10 plotEmpirical.py