import numpy as np
from hmm_utils import backward_algorithm, _nstep_log_trans_prob
from scipy.special import logsumexp
import argparse
import os
from scipy.stats import norm, beta

def parse_args():
	"""Define the Arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--times',type=str, default=None)
	parser.add_argument('--popFreq',type=float,default=None)

	parser.add_argument('--ancientSamps',type=str,default=None)
	parser.add_argument('--ancientHaps',type=str,default=None)
	parser.add_argument('--out',type=str,default=None)

	parser.add_argument('--N',type=float,default=None)
	parser.add_argument('--coal',type=str,default=None,help='path to Relate .coal file. Negates --N option.')

	parser.add_argument('--tCutoff',type=float,default=None)
	parser.add_argument('--timeBins',type=str,default=None, nargs='+')
	parser.add_argument('--sMax',type=float,default=0.1)
	parser.add_argument('--df',type=int,default=450)
	parser.add_argument('--noAlleleTraj', default=False, action='store_true', help='whether to compute the posterior allele frequency trajectory or not.')
	parser.add_argument('--integration_points', type=int, default = -1)
	parser.add_argument('--h', type=float, default = 0.5)

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
	if "," not in Lines[0]:
		DerivedSampleTimes = Lines[0]
		if DerivedSampleTimes == "\n":
			derstimes = np.array([])
		else:
			DerivedSampleTimes = DerivedSampleTimes.split(";")
			derstimes = -1.0 * np.ones(len(DerivedSampleTimes))
			for i in range(len(derstimes)):
				derstimes[i] = float(DerivedSampleTimes[i])
		AncestralSampleTimes = Lines[1]
		if AncestralSampleTimes == "\n":
			ancstimes = np.array([])
		else:
			AncestralSampleTimes = AncestralSampleTimes.split(";")
			ancstimes = -1.0 * np.ones(len(AncestralSampleTimes))
			
			for i in range(len(ancstimes)):
				ancstimes[i] = float(AncestralSampleTimes[i])
		Lines = Lines[2:]
	else:
		derstimes = np.array([])
		ancstimes = np.array([])
	
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
	return(np.array([dertimes,anctimes]),derstimes,ancstimes)
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
		times,derSampledTimes, ancSampledTimes = load_times(args.times)
	else:
		times,derSampledTimes, ancSampledTimes = np.zeros((2,0,0)),np.array([]),np.array([])
	currFreq = args.popFreq

	# load ancient samples/genotype likelihoods
	if args.ancientSamps != None:
		ancientGLs = np.genfromtxt(args.ancientSamps,delimiter=' ')
	else:
		ancientGLs = np.zeros((0,4))

	# load ancient haploid genotype likelihoods
	if args.ancientHaps != None:
		ancientHapGLs = np.genfromtxt(args.ancientHaps,delimiter=' ')
	else:
		ancientHapGLs = np.zeros((0,3))

	tCutoff = args.tCutoff

	epochs = np.arange(0.0,tCutoff,int(1))
	# loading population size trajectory
	if args.coal != None:
		Nepochs = np.genfromtxt(args.coal,skip_header=1,skip_footer=1)
		N = 1/np.genfromtxt(args.coal,skip_header=2)[2:-1] # should be 1
		N = np.array(list(N)+[N[-1]])
		Ne = N[np.digitize(epochs,Nepochs)-1]
	else:
		if args.N == None:
			print("Must specify either N or .coal file")
			exit()
		Ne = args.N * np.ones(int(tCutoff))
	# set up freq bins
	freqs = beta.ppf(np.linspace(0.0, 1.0, args.df), 0.5, 0.5)
	freqs[0] = 1e-12
	freqs[len(freqs) - 1] = 1 - 1e-12
	logfreqs = np.log(freqs)
	log1minusfreqs = np.log(np.subtract(1,freqs))
	freqs[0] = 0.0
	freqs[len(freqs) - 1] = 1.0

	# load time bins (for defining selection epochs)
	if args.timeBins != None:
		breakppoint = []
		for i in range(len(args.timeBins)):
			breakppoint.append(float(round(float(args.timeBins[i]))))
		timeBins = np.array([0.0] + breakppoint + [tCutoff + 20.0])
	else:
		timeBins = np.array([0.0,tCutoff])
	return timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,logfreqs,log1minusfreqs,derSampledTimes,ancSampledTimes,args.h

def likelihood_wrapper(theta,timeBins,N,h,freqs,times, logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,Weights = []):
	S = theta
	Sprime = np.concatenate((S,[0.0]))
	if np.any(np.abs(Sprime) > sMax):
		return 1e+100 *( 10**(np.max(np.abs(Sprime)))/sMax)
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
		tranmatrix = _nstep_log_trans_prob(N[0],sel[0],freqs,z_bins,z_logcdf,z_logsf,h) # this only handles the precompute case, just use initial values
		if len(np.unique(N)) + len(np.unique(sel)) == 2:
			precompute = 1
		for i in range(M):
			betaMat = backward_algorithm(sel,times[:,:,i],derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
			#np.savetxt(str(i)+"_beta.txt", betaMat, delimiter=",")
			logl = logsumexp(betaMat[-2,:])
			loglrs[i] = logl-Weights[i]
		logl = -1 * (-np.log(M) + logsumexp(loglrs))
	else:
		betaMat = backward_algorithm(sel,t,derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,np.zeros((len(freqs),len(freqs))),noCoals=noCoals,precomputematrixboolean=0,currFreq=currFreq)
		logl = -logsumexp(betaMat[-2,:])
	strr = ""
	for i in range(len(theta)):
		strr = strr + str(theta[i]) + ","
	functionvals.write(strr + str(-logl) + "\n")
	return logl

if __name__ == "__main__":
	args = parse_args()
	if args.times == None and args.ancientSamps == None and args.ancientHaps == None:
		print('You need to supply coalescence times (--times) and/or ancient samples (--ancientSamps)')
	
	# load data and set up model
	#should delete the previous verion
	if os.path.exists(args.out+"_tempfile.txt"):
		os.remove(args.out+"_tempfile.txt")
	functionvals = open(args.out+"_tempfile.txt", "a")
	sMax = args.sMax
	timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,logfreqs,log1minusfreqs,derSampledTimes,ancSampledTimes,h = load_data(args)
	# read in global Phi(z) lookups
	linearspacing = np.linspace(0.0, 1.0, 2000)
	linearspacing[0] = 1e-10
	linearspacing[len(linearspacing) - 1] = 1 - 1e-10

	z_bins = norm.ppf(linearspacing)
	z_logcdf = norm.cdf(z_bins)
	z_logsf = z_logcdf # can delete this later.

	Ne *= 1/2
	noCoals = int(noCoals)

	# optimize over selection parameters
	T = len(timeBins)
	S0 = 0.0 * np.ones(T-1)

	M = times.shape[2]
	Weights = np.zeros(M)
	precompute = 0
	tranmatrix = _nstep_log_trans_prob(Ne[0],0.0,freqs,z_bins,z_logcdf,z_logsf,h) # this only handles the precompute case, just use initial values
	if len(np.unique(Ne)) == 1:
		precompute = 1
	for i in range(M):
		betaMatl0 = backward_algorithm(np.zeros(len(Ne)),times[:,:,i],derSampledTimes,ancSampledTimes,epochs,Ne,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
		Weights[i] = logsumexp(betaMatl0[-2,:])

	minargs = (timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,Weights)

	for iiiiii in np.arange(0.014,0.0400001, 0.0004):
			print(iiiiii, -likelihood_wrapper([iiiiii],timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,Weights))