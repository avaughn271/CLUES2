from matplotlib import pyplot as plt
import matplotlib.cm as cmm

import numpy as np
import argparse
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser()
args = parser.parse_args()
ColorIntensity = []
for i in range(3):
    ColorIntensity.append(1 - float(i)/3)
ConfidenceIntervals = [0.5,0.75,0.95]

freqs =  np.loadtxt("output_freqs.txt", delimiter=",", dtype=float)
logpost = np.loadtxt("output_post.txt", delimiter=",", dtype=float)
epochs = np.linspace(0,logpost.shape[1],logpost.shape[1] + 1)
f,ax = plt.subplots(1,1)
f.set_size_inches(20,10)
EXPPOST = logpost # rows are frequency, columns are time interbals

MATRIXTOPLOT = np.zeros((logpost.shape[0],logpost.shape[1]))

for indexx in range(len(ConfidenceIntervals)):
    lower =  np.loadtxt("Lower" + str(ConfidenceIntervals[indexx])+ ".txt", delimiter=",", dtype=float)
    upper =  np.loadtxt("Upper" + str(ConfidenceIntervals[indexx])+ ".txt", delimiter=",", dtype=float)
    print(lower)
    for timeinterval in epochs:
        LOWERMAX = 10.0
        UPPERMAX = 10.0
        currentlow = -1
        currenthigh = -1
        TimeSLICE = EXPPOST[:,int(timeinterval)]
        for frequency in range(0, len(freqs)):
            fff = freqs[frequency]
            if abs(fff - lower[int(timeinterval)]) < LOWERMAX:
                 LOWERMAX = abs(fff - lower[int(timeinterval)])
                 currentlow = frequency
            if abs(fff - upper[int(timeinterval)]) < UPPERMAX:
                 UPPERMAX = abs(fff - upper[int(timeinterval)])
                 currenthigh = frequency

        for j in range(MATRIXTOPLOT.shape[0]):
            if j >= currentlow and j <= currenthigh:
                MATRIXTOPLOT[j, int(timeinterval)] = np.max([MATRIXTOPLOT[j, int(timeinterval)], ColorIntensity[indexx]])
        if timeinterval == epochs[len(epochs) - 2]: break

plt.pcolormesh(epochs[:-1],freqs, MATRIXTOPLOT)
plt.axis((0,len(epochs[:-1]),0,1.0))
plt.ylabel('Derived allele frequency',fontsize=20)
plt.xlabel('Generations before present',fontsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

cmap = cmm.get_cmap()
handless = []
for i in range(len(ConfidenceIntervals)):
    handless.append(mpatches.Patch(color=cmap(ColorIntensity[i]), label= str(ConfidenceIntervals[i] * 100) + '% Posterior Interval'))

ax.legend(handles=handless, loc='upper right')

plt.savefig('empirical.png',format='png', dpi=800, bbox_inches='tight')
