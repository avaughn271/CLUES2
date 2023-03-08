from matplotlib import pyplot as plt
import matplotlib

import numpy as np
import argparse
from matplotlib import collections  as mc

parser = argparse.ArgumentParser()
parser.add_argument('--freqs',type=str)
parser.add_argument('--post',type=str)
parser.add_argument('--figure',type=str)
args = parser.parse_args()

freqs =  np.loadtxt(args.freqs, delimiter=",", dtype=float)
logpost = np.loadtxt(args.post, delimiter=",", dtype=float)
epochs = np.linspace(0,logpost.shape[1],logpost.shape[1] + 1)
f,ax = plt.subplots(1,1)
f.set_size_inches(20,10)
EXPPOST = np.exp(logpost) # rows are frequency, columns are time interbals
Max = (np.max(EXPPOST[:,logpost.shape[1]-1])) * 5
norm = matplotlib.colors.Normalize(vmin=0.0, vmax=Max)
plt.pcolormesh(epochs[:-1],freqs,EXPPOST[:,:],norm=norm)
plt.axis((0,len(epochs[:-1]),0,1.0))
plt.ylabel('Derived allele frequency',fontsize=20)
plt.xlabel('Generations before present',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

cbar = plt.colorbar(norm=norm)
cbar.ax.set_ylabel('Posterior probability\n\n',rotation=270,fontsize=20,labelpad=40)
cbar.ax.tick_params(labelsize=18)

LINES = []

for timeinterval in epochs:
    if timeinterval == epochs[len(epochs)-1]: break
    tentativespan = 10000000.0
    currentlow = -1
    currenthigh = -1
    TimeSLICE = EXPPOST[:,int(timeinterval)]
    for lowerfrequencyindex in range(len(freqs)):
        for higherfrequencyindex in range(lowerfrequencyindex + 1, len(freqs)):
            possiblesum = np.sum(TimeSLICE[lowerfrequencyindex:higherfrequencyindex])
            rangeofsum = freqs[higherfrequencyindex] - freqs[lowerfrequencyindex]
            if possiblesum >= 0.95 and rangeofsum < tentativespan:
                tentativespan = rangeofsum
                currentlow = freqs[lowerfrequencyindex]
                currenthigh = freqs[higherfrequencyindex]
            elif rangeofsum > tentativespan:
                break
    LINES.append([(timeinterval,currentlow),(timeinterval+1,currentlow)])
    LINES.append([(timeinterval,currenthigh) ,(timeinterval+1,currenthigh)])

#iterate through interval starts and interval ends, find the smallest interval that contains
#at least 0.95 of the probability.
lc = mc.LineCollection(LINES, colors=(1, 0, 0, 1), linewidths=2)
ax.add_collection(lc)


plt.savefig('%s.%s'%(args.figure,'png'),format='png', dpi=800, bbox_inches='tight')
