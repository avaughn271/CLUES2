

python3 inference.py --popFreq 0.55  --ancientSamps Samples1.txt --N 30000 --out outt --tCutoff 600 --df 300

head -n 1 outt_post.txt

head outt_inference.txt

python3 plot_traj.py --freqs outt_freqs.txt --post outt_post.txt --figure figoutput

