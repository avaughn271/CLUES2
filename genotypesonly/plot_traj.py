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

plt.pcolormesh(epochs[:-1],freqs,np.exp(logpost)[:,:])
plt.axis((0,len(epochs[:-1]),0,1.0))
plt.ylabel('Allele frequency',fontsize=20)
plt.xlabel('Generations before present',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

cbar = plt.colorbar()
cbar.ax.set_ylabel('Posterior probability\n\n',rotation=270,fontsize=20,labelpad=40)
cbar.ax.tick_params(labelsize=18)

plt.savefig('%s.%s'%(args.figure,'png'),format='png', dpi=800, bbox_inches='tight')
