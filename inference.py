import numpy as np
from hmm_utils import backward_algorithm
from scipy.special import logsumexp
from scipy.optimize import minimize_scalar
import argparse
import os

def parse_args():
	"""Define the Arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--times',type=str, default=None)
	parser.add_argument('--popFreq',type=float,default=None)

	parser.add_argument('--ancientSamps',type=str,default=None)
	parser.add_argument('--out',type=str,default=None)

	parser.add_argument('-N','--N',type=float,default=10**4)
	parser.add_argument('-coal','--coal',type=str,default=None,help='path to Relate .coal file. Negates --N option.')

	parser.add_argument('--tCutoff',type=float,default=1000)
	parser.add_argument('--timeBins',type=str,default=None)
	parser.add_argument('--sMax',type=float,default=0.1)
	parser.add_argument('--df',type=int,default=450)
	parser.add_argument('--noAlleleTraj', default=False, action='store_true', help='whether to compute the posterior allele frequency trajectory or not.')

	return parser.parse_args()

def load_times(readtimes):
	file1 = open(readtimes, 'r')
	Lines = file1.readlines()
	
	M = int(len(Lines) / 2)
	ntot = len(Lines[0].split(",")) + len(Lines[1].split(","))
	dertimes = -1.0 * np.ones((ntot,M))
	anctimes = -1.0 * np.ones((ntot,M))

	for m in range(M):
		der = Lines[2 * m].split(",")
		anc = Lines[2 * m + 1].split(",")
		for i in range(len(der)):
			dertimes[i,m] = float(der[i])
		for i in range(len(anc)):
			anctimes[i,m] = float(anc[i])
	return(np.array([dertimes,anctimes]),np.array([]),np.array([]))

def load_data(args):
	noCoals = (args.times == None)
	times,derSampledTimes, ancSampledTimes = load_times(args.times)
	currFreq = args.popFreq

	ancientGLs = np.zeros((0,4))

	# load ancient haploid genotype likelihoods
	ancientHapGLs = np.zeros((0,3))

	tCutoff = args.tCutoff

	epochs = np.arange(0.0,tCutoff,int(1))

	Ne = args.N * np.ones(int(tCutoff))
	# set up freq bins
	beta05quantiles = np.genfromtxt(os.path.dirname(__file__) + '/utils/Beta05Quantiles.txt')
	freqs = np.quantile(beta05quantiles, np.linspace(0.0, 1.0, args.df))
	logfreqs = np.log(freqs)
	log1minusfreqs = np.log(np.subtract(1,freqs))

	timeBins = np.array([0.0,tCutoff])
	return timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,logfreqs,log1minusfreqs,derSampledTimes,ancSampledTimes

def likelihood_wrapper(theta,timeBins,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes, Weights = []):
	S = theta
	Sprime = np.concatenate((S,[0.0]))
	if np.any(np.abs(Sprime) > sMax):
		return 1e+100 *( 10**(np.max(np.abs(Sprime)))/sMax)
	sel = Sprime[np.digitize(epochs,timeBins,right=False)-1]
	t = times[:,:,0]

	betaMat = backward_algorithm(sel,t,derSampledTimes,ancSampledTimes,epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,np.zeros((len(freqs),len(freqs))),noCoals=noCoals,precomputematrixboolean=0,currFreq=currFreq)
	logl = -logsumexp(betaMat[-2,:])
	return logl

def likelihood_wrapper_scalar(theta,timeBins,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,Weights = []):
	return(likelihood_wrapper([theta - 1.0],timeBins,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,Weights))
#added +1 and -1 in order to get better convergence properties.

if __name__ == "__main__":
	args = parse_args()

	sMax = args.sMax
	timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,logfreqs,log1minusfreqs,derSampledTimes,ancSampledTimes = load_data(args)
	# read in global Phi(z) lookups
	z_bins = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_bins.txt')
	z_logcdf = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_logcdf.txt')
	z_logsf = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_logsf.txt')

	Ne *= 1/2
	noCoals = int(noCoals)
	rangerr = 8
	K = [x/(100 * rangerr) for x in range(98*rangerr, 102*rangerr)]
	print(K)
	for k in K:
		print(k-1.0, likelihood_wrapper_scalar(k,timeBins,Ne,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes))