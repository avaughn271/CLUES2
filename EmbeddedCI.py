from matplotlib import pyplot as plt

import numpy as np
import argparse

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

MATRIXTOPLOT = np.zeros((EXPPOST.shape[0],EXPPOST.shape[1]))


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
            if possiblesum >= 0.5 and rangeofsum < tentativespan:
                tentativespan = rangeofsum
                currentlow = lowerfrequencyindex
                currenthigh = higherfrequencyindex
                #if lowerfrequencyindex == higherfrequencyindex: print(TimeSLICE)
            elif rangeofsum > tentativespan:
                break
    for j in range(MATRIXTOPLOT.shape[0]):
        if j >= currentlow and j <= currenthigh:
            MATRIXTOPLOT[j, int(timeinterval)] = np.max([MATRIXTOPLOT[j, int(timeinterval)], 1.0])


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
            if possiblesum >= 0.75 and rangeofsum < tentativespan:
                tentativespan = rangeofsum
                currentlow = lowerfrequencyindex
                currenthigh = higherfrequencyindex
                #if lowerfrequencyindex == higherfrequencyindex: print(TimeSLICE)
            elif rangeofsum > tentativespan:
                break
    for j in range(MATRIXTOPLOT.shape[0]):
        if j >= currentlow and j <= currenthigh:
            MATRIXTOPLOT[j, int(timeinterval)] = np.max([MATRIXTOPLOT[j, int(timeinterval)], 0.75])




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
                currentlow = lowerfrequencyindex
                currenthigh = higherfrequencyindex
                #if lowerfrequencyindex == higherfrequencyindex: print(TimeSLICE)
            elif rangeofsum > tentativespan:
                break
    for j in range(MATRIXTOPLOT.shape[0]):
        if j >= currentlow and j <= currenthigh:
            MATRIXTOPLOT[j, int(timeinterval)] = np.max([MATRIXTOPLOT[j, int(timeinterval)], 0.5])






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
            if possiblesum >= 0.99 and rangeofsum < tentativespan:
                tentativespan = rangeofsum
                currentlow = lowerfrequencyindex
                currenthigh = higherfrequencyindex
                #if lowerfrequencyindex == higherfrequencyindex: print(TimeSLICE)
            elif rangeofsum > tentativespan:
                break
    for j in range(MATRIXTOPLOT.shape[0]):
        if j >= currentlow and j <= currenthigh:
            MATRIXTOPLOT[j, int(timeinterval)] = np.max([MATRIXTOPLOT[j, int(timeinterval)], 0.25])





plt.pcolormesh(epochs[:-1],freqs, MATRIXTOPLOT)
plt.axis((0,len(epochs[:-1]),0,1.0))
plt.ylabel('Derived allele frequency',fontsize=20)
plt.xlabel('Generations before present',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

cbar = plt.colorbar()
cbar.ax.set_ylabel('Posterior probability\n\n',rotation=270,fontsize=20,labelpad=40)
cbar.ax.tick_params(labelsize=18)

#iterate through interval starts and interval ends, find the smallest interval that contains
#at least 0.95 of the probability.

plt.savefig('%s.%s'%(args.figure,'png'),format='png', dpi=800, bbox_inches='tight')
