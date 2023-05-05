import numpy as np
from hmm_utils import backward_algorithm, forward_algorithm, _nstep_log_trans_prob
from scipy.special import logsumexp
from scipy.optimize import minimize, minimize_scalar
import argparse
import os
from scipy.stats import chi2

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
	"""Load in the coalescent times
    INPUT: readtimes - this is the file name of coalescence times
           times are given on an absolute scale.
	   Line 1 is coalescence times of derived lineages including the mixed lineage!!!!
	   Line 2 is coalescence times of ancestral lineages
	   This is repeated for each important sample
	   If there are A derived coalescences, there are A derived samples
	   If there are B ancestral coalescences, there are B + 1 ancestral samples

    OUTPUT: A numpy array of dimension (2,A+B, numberimportancesamples)
    #NEED TO DOUBLE CHECK WITH PAPER FOR HOW MIXED ANCESTRAL IS TAKEN INTO ACCOUNT AND HOW
    YOU FILTER TIMES.
    For the first dimension of the output, you populate it with the first B of the B coalescences
    For the second dimension of the output, you populate it with the first A-1 of the A coalescences
    All of the other entries are -1.
    """
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

	return(np.array([dertimes,anctimes]))
	#locus time is an array that has dimensions 2 by (total number of leaves) by (number of importance samples)
	#The first row corresponds to the derived alleles. The columns are populated by daf-1 and n-daf-1 entries each, and are then -1 below this value.

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
	noCoals = (args.times == None)
	if not noCoals:
		times = load_times(args.times)
	else:
		times = np.zeros((2,0,0))
	currFreq = args.popFreq

	# load ancient samples/genotype likelihoods
	if args.ancientSamps != None:
		ancientGLs = np.genfromtxt(args.ancientSamps,delimiter=' ')
	else:
		ancientGLs = np.zeros((0,4))

	# load ancient haploid genotype likelihoods
	ancientHapGLs = np.zeros((0,3))

	tCutoff = args.tCutoff

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
	logfreqs = np.log(freqs)
	log1minusfreqs = np.log(np.subtract(1,freqs))

	# load time bins (for defining selection epochs)
	if args.timeBins != None:
		timeBins = np.genfromtxt(args.timeBins)
	else:
		timeBins = np.array([0.0,tCutoff])
	return timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,logfreqs,log1minusfreqs

def likelihood_wrapper(theta,timeBins,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax, Weights = []):
	S = theta
	Sprime = np.concatenate((S,[0.0]))
	if np.any(np.abs(Sprime) > sMax):
		return np.inf
	sel = Sprime[np.digitize(epochs,timeBins,right=False)-1]
	tShape = times.shape
	if tShape[2] == 0:
		t = np.zeros((2,0))
		importanceSampling = False
	elif tShape[2] == 1:
		t = times[:,:,0]
		importanceSampling = False
	else:
		importanceSampling = True
	if importanceSampling:
		M = tShape[2]
		loglrs = np.zeros(M)
		precompute = 0
		tranmatrix = _nstep_log_trans_prob(N[0],sel[0],freqs,z_bins,z_logcdf,z_logsf) # this only handles the precompute case, just use initial values
		if len(np.unique(N)) + len(np.unique(sel)) == 2:
			precompute = 1
		for i in range(M):
			betaMat = backward_algorithm(sel,times[:,:,i],epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
			#np.savetxt(str(i)+"_beta.txt", betaMat, delimiter=",")
			logl = logsumexp(betaMat[-2,:])
			loglrs[i] = logl-Weights[i]
		logl = -1 * (-np.log(M) + logsumexp(loglrs))
	else:
		betaMat = backward_algorithm(sel,t,epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,np.zeros((len(freqs),len(freqs))),noCoals=noCoals,precomputematrixboolean=0,currFreq=currFreq)
		logl = -logsumexp(betaMat[-2,:])
	return logl

def likelihood_wrapper_scalar(theta,timeBins,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax, Weights = []):
	return(likelihood_wrapper([theta - 1.0],timeBins,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax, Weights))
#added +1 and -1 in order to get better convergence properties.

def traj_wrapper(theta,timeBins,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax, Weights = []):
	S = theta
	Sprime = np.concatenate((S,[0.0]))
	if np.any(np.abs(Sprime) > sMax):
		print('WARNING: selection coefficient exceeds bounds. Maybe change --sMax?')
		return np.inf

	sel = Sprime[np.digitize(epochs,timeBins,right=False)-1]
	T = len(epochs)
	F = len(freqs)
	tShape = times.shape
	if tShape[2] == 0:
		t = np.zeros((2,0))
		importanceSampling = False
	elif tShape[2] == 1:
		t = times[:,:,0]
		importanceSampling = False
	else:
		importanceSampling = True

	if importanceSampling:
		M = tShape[2]
		loglrs = np.zeros(M)
		postBySamples = np.zeros((F,T-1,M))
		tranmatrix = _nstep_log_trans_prob(N[0],sel[0],freqs,z_bins,z_logcdf,z_logsf) # this only handles the precompute case, just use initial values
		precompute = 0
		if len(np.unique(N)) + len(np.unique(sel)) == 2:
			precompute = 1
		for i in range(M):
			betaMat = backward_algorithm(sel,times[:,:,i],epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
			alphaMat = forward_algorithm(sel,times[:,:,i],epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals)
			logl = logsumexp(betaMat[-2,:])
			loglrs[i] = logl - Weights[i]
			postBySamples[:,:,i] = (alphaMat[1:,:] + betaMat[:-1,:]).transpose()
		post = logsumexp(loglrs + postBySamples,axis=2)
		post -= logsumexp(post,axis=0)

	else:
		post = np.zeros((F,T))
		betaMat = backward_algorithm(sel,t,epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,np.zeros((F,F)),noCoals=noCoals,precomputematrixboolean=0,currFreq=currFreq)
		alphaMat = forward_algorithm(sel,t,epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals)
		post = (alphaMat[1:,:] + betaMat[:-1,:]).transpose()
		post -= logsumexp(post,axis=0)
	return post

if __name__ == "__main__":
	args = parse_args()
	if args.times == None and args.ancientSamps == None:
		print('You need to supply coalescence times (--times) and/or ancient samples (--ancientSamps)')
	
	# load data and set up model
	sMax = args.sMax	
	timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,logfreqs,log1minusfreqs = load_data(args)
	# read in global Phi(z) lookups
	z_bins = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_bins.txt')
	z_logcdf = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_logcdf.txt')
	z_logsf = np.genfromtxt(os.path.dirname(__file__) + '/utils/z_logsf.txt')

	Ne *= 1/2
	noCoals = int(noCoals)

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
	ImpSamp = False
	if times.shape[2] > 1:
		print('\t(Importance sampling with M = %d samples)'%(times.shape[2]))
		print()
		ImpSamp = True
	if not ImpSamp: # to account for whether we return the likelihood or the log likelihood

		logL0 = likelihood_wrapper(S0,timeBins,Ne,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax)

		minargs = (timeBins,Ne,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax)

		if len(S0) == 1:
			res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.9,1.0,1.1],args=minargs, method = "Brent", tol = 1e-4))
			S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
			L = res.fun
			print(S,L)
		else:
			res = minimize(likelihood_wrapper, S0, args=minargs, options=opts, method='Nelder-Mead')
			S = res.x
			L = res.fun

		toprint = '%.4f'%(-L+logL0)
		numericloglik = -L+logL0

		Weights = []
	else:
		M = times.shape[2]
		Weights = np.zeros(M)
		precompute = 0
		tranmatrix = _nstep_log_trans_prob(Ne[0],0.0,freqs,z_bins,z_logcdf,z_logsf) # this only handles the precompute case, just use initial values
		if len(np.unique(Ne)) == 1:
			precompute = 1
		for i in range(M):
			betaMatl0 = backward_algorithm(np.zeros(len(Ne)),times[:,:,i],epochs,Ne,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
			Weights[i] = logsumexp(betaMatl0[-2,:])

		minargs = (timeBins,Ne,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax, Weights)

		if len(S0) == 1:
			res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.9,1.0,1.1],args=minargs, method = "Brent", tol = 1e-4))
			S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
			L = res.fun
			print(S,L)
		else:
			res = minimize(likelihood_wrapper, S0, args=minargs, options=opts, method='Nelder-Mead')
			S = res.x
			L = res.fun
		numericloglik = -L
		toprint = '%.4f'%(-L)

	FirstLine = "logLR" + "\t" + "log10(p-value)"
	epochnum = 1
	degreesoffreedom = len(timeBins) - 1
	toprint = toprint + "\t" + '%.2f'%((chi2.logsf(numericloglik + numericloglik, degreesoffreedom ) ) / np.log(10) )
	for s,t,u in zip(S,timeBins[:-1],timeBins[1:]):
		toprint = toprint + "\t" + '%d'%(t)
		toprint = toprint + "\t" + '%d'%(u)
		toprint = toprint + "\t" + '%.5f'%(s)
		FirstLine = FirstLine + "\t" + "Epoch" + str(epochnum) + "_start" + "\t" +  "Epoch" + str(epochnum)   + "_end"  + "\t" +  "SelectionMLE" + str(epochnum)
		epochnum = epochnum + 1
	f = open(args.out+"_inference.txt", "w+")
	f.writelines(FirstLine  + "\n" + toprint + "\n")
	f.close()

	if not args.noAlleleTraj:
		# infer trajectory @ MLE of selection parameter
		post = traj_wrapper(S,timeBins,Ne,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax, Weights)
		np.savetxt(args.out+"_freqs.txt", freqs, delimiter=",")
		np.savetxt(args.out+"_post.txt", post, delimiter=",")