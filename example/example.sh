#This is the basic command to run CLUES2 on samples of ARGs
python PATH/CLUES2/inference.py --popFreq  0.8 --times InputTimes.txt --N 30000 --out output1 --tCutoff 1000

#This is the basic command to run CLUES2 on ancient genotype data
python PATH/CLUES2/inference.py --popFreq  0.8 --ancientSamps exampleAncientSamps.txt --N 30000 --out output2 --tCutoff 1000

#We can also plot the reconstructed derived allele frequencies of the above commands
python PATH/CLUES2/plot_traj.py --freqs  output1_freqs.txt --post output1_post.txt --figure ExampleFigure1  --generation_time 28.0
python PATH/CLUES2/plot_traj.py --freqs  output2_freqs.txt --post output2_post.txt --figure ExampleFigure2  --generation_time 28.0