import numpy as np
from hmm_utils import forward_algorithm
from hmm_utils import backward_algorithm
from scipy.special import logsumexp
from scipy.optimize import minimize
import argparse
import os

def parse_args():
	"""Define the Arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--popFreq',type=float,default=None)

	parser.add_argument('--ancientSamps',type=str,default=None)
	parser.add_argument('--out',type=str,default=None)

	parser.add_argument('-N','--N',type=float,default=10**4)
	parser.add_argument('-coal','--coal',type=str,default=None,help='path to Relate .coal file. Negates --N option.')

	parser.add_argument('--tCutoff',type=float,default=1000)
	parser.add_argument('--timeBins',type=str,default=None)
	parser.add_argument('--sMax',type=float,default=0.1)
	parser.add_argument('--df',type=int,default=400)
	return parser.parse_args()

def load_data(args):
	"""Takes in arguments.
	Returns:
	timeBins - if tcutoff given, this is [0.0, tcutoff]. If bins specified, use them instead.
	times - the output of load_times. or 0 if none specified
	epochs - [0, ..., tcutoff - 1]
	Ne - [Ne] * tcutoff or equivalent if coal is given. This is haploid size (larger one)!
	freqs - array of increasing numbers from 0 to 1 exclusive. Length is args.df 
	ancientGLs - matrix of log proabilities of observations. First column is times. 
	             Second is homo anc. Last is homo derived. All 0 if not specified
	ancientHapGLs - matrix of zeros.
	noCoals - True if times not specified. False otherwise.
	currFreq - user specified modern derived allele frequency.
	"""
	currFreq = args.popFreq

	# load ancient samples/genotype likelihoods
	if args.ancientSamps != None:
		ancientGLs = np.genfromtxt(args.ancientSamps,delimiter=' ')
	else:
		ancientGLs = np.zeros((0,4))

	# load ancient haploid genotype likelihoods
	ancientHapGLs = np.zeros((0,3))

	tCutoff = np.max(ancientGLs[:,0])+1.0

	epochs = np.arange(0.0,tCutoff,int(1))
	# loading population size trajectory
	if args.coal != None:
		Nepochs = np.genfromtxt(args.coal,skip_header=1,skip_footer=1)
		N = 0.5/np.genfromtxt(args.coal,skip_header=2)[2:-1]
		N = np.array(list(N)+[N[-1]])
		Ne = N[np.digitize(epochs,Nepochs)-1]
	else:
		Ne = args.N * np.ones(int(tCutoff))
	# set up freq bins
	beta05quantiles = np.genfromtxt(os.path.dirname(__file__) + '/utils/Beta05Quantiles.txt')
	freqs = np.quantile(beta05quantiles, np.linspace(0.0, 1.0, args.df))
	# load time bins (for defining selection epochs)
	if args.timeBins != None:
		timeBins = np.genfromtxt(args.timeBins)
	else:
		timeBins = np.array([0.0,tCutoff])
	return timeBins,epochs,Ne,freqs,ancientGLs,ancientHapGLs,currFreq

def likelihood_wrapper(theta,timeBins,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,currFreq,sMax):
	Sprime = np.concatenate((theta,[0.0]))
	if np.any(np.abs(Sprime) > sMax):
		return np.inf
	sel = Sprime[np.digitize(epochs,timeBins,right=False)-1]

	betaMat = backward_algorithm(sel,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,currFreq=currFreq)
	logl = -logsumexp(betaMat[-1,:])
	return logl

def traj_wrapper(theta,timeBins,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,currFreq,sMax):
	S = theta
	Sprime = np.concatenate((S,[0.0]))
	if np.any(np.abs(Sprime) > sMax):
		print('WARNING: selection coefficient exceeds bounds. Maybe change --sMax?')
		return np.inf

	sel = Sprime[np.digitize(epochs,timeBins,right=False)-1]

	betaMat = backward_algorithm(sel,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,currFreq=currFreq)
	alphaMat = forward_algorithm(sel,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs)
	post = (alphaMat + betaMat).transpose()
	post -= logsumexp(post,axis=0)
	return post

if __name__ == "__main__":
	args = parse_args()
	if args.ancientSamps == None:
		print('You need to supply ancient samples (--ancientSamps)')
	# load data and set up model
	sMax = args.sMax
	timeBins,epochs,Ne,freqs,ancientGLs,ancientHapGLs,currFreq = load_data(args)
	# read in global Phi(z) lookups
	z_bins = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_bins.txt')
	z_logcdf = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_logcdf.txt')
	z_logsf = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_logsf.txt')

	Ne *= 1/2

	# optimize over selection parameters
	T = len(timeBins)
	S0 = 0.0 * np.ones(T-1)
	opts = {'xatol':1e-4}

	if T == 2:
		Simplex = np.reshape(np.array([-0.05,0.05]),(2,1))
	elif T > 2:
		Simplex = np.zeros((T,T-1))
		for i in range(Simplex.shape[1]):
			Simplex[i,:] = -0.01
			Simplex[i,i] = 0.01
		Simplex[-1,:] = 0.01
	else:
		raise ValueError

	opts['initial_simplex']=Simplex
	logL0 = likelihood_wrapper(S0,timeBins,Ne,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,currFreq,sMax)

	minargs = (timeBins,Ne,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,currFreq,sMax)
	res = minimize(likelihood_wrapper, S0, args=minargs, options=opts, method='Nelder-Mead')

	toprint = []

	toprint.append('logLR: %.4f'%(-res.fun+logL0) + "\n")
	toprint.append('Epoch\tSelection MLE'+ "\n")
	for s,t,u in zip(res.x,timeBins[:-1],timeBins[1:]):
		toprint.append('%d-%d\t%.5f'%(t,u,s)+ "\n")

	# infer trajectory @ MLE of selection parameter
	post = traj_wrapper(res.x,timeBins,Ne,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,currFreq,sMax)
	
	f = open(args.out+"_inference.txt", "w+")
	f.writelines(toprint)
	f.close()
	np.savetxt(args.out+"_post.txt", post, delimiter=",") #print(i,np.sum(freqs * np.exp(post[:,i])))
	np.savetxt(args.out+"_freqs.txt", freqs, delimiter=",") #print(i,np.sum(freqs * np.exp(post[:,i])))
