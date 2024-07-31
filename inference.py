import numpy as np
from hmm_utils import backward_algorithm, forward_algorithm, _nstep_log_trans_prob
from scipy.special import logsumexp
from scipy.optimize import minimize, minimize_scalar
import argparse
import os
from scipy.stats import chi2, norm, beta
from numpy.random import normal
from scipy.stats import multivariate_normal

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
		for i in range(len(ancientGLs)):
			liksum = ancientGLs[i,1] + ancientGLs[i,2] + ancientGLs[i,3]
			if liksum < 1.001 and liksum > 0.999:
				ancientGLs[i, 1:4] = np.log(ancientGLs[i, 1:4] +  np.array([1e-45,1e-45,1e-45]) )
	else:
		ancientGLs = np.zeros((0,4))

	# load ancient haploid genotype likelihoods
	if args.ancientHaps != None:
		ancientHapGLs = np.genfromtxt(args.ancientHaps,delimiter=' ')
		for i in range(len(ancientHapGLs)):
			liksum = ancientHapGLs[i,1] + ancientHapGLs[i,2]
			if liksum < 1.001 and liksum > 0.999:
				ancientHapGLs[i, 1:3] = np.log(ancientHapGLs[i, 1:3] +  np.array([1e-45,1e-45]) )
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

def likelihood_wrapper_scalar(theta,timeBins,N,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,Weights = []):
	return(likelihood_wrapper([theta - 1.0],timeBins,N,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,Weights))
#added +1 and -1 in order to get better convergence properties.

def traj_wrapper(theta,timeBins,N,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,Weights = []):
	S = theta
	Sprime = np.concatenate((S,[0.0]))

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
		tranmatrix = _nstep_log_trans_prob(N[0],sel[0],freqs,z_bins,z_logcdf,z_logsf,h) # this only handles the precompute case, just use initial values
		precompute = 0
		if len(np.unique(N)) + len(np.unique(sel)) == 2:
			precompute = 1
		for i in range(M):
			betaMat = backward_algorithm(sel,times[:,:,i],derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
			alphaMat = forward_algorithm(sel,times[:,:,i],derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals)
			logl = logsumexp(betaMat[-2,:])
			loglrs[i] = logl - Weights[i]
			postBySamples[:,:,i] = (alphaMat[1:,:] + betaMat[:-1,:]).transpose()
		post = logsumexp(loglrs + postBySamples,axis=2)
		post -= logsumexp(post,axis=0)

	else:
		post = np.zeros((F,T))
		betaMat = backward_algorithm(sel,t,derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,np.zeros((F,F)),noCoals=noCoals,precomputematrixboolean=0,currFreq=currFreq)
		alphaMat = forward_algorithm(sel,t,derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals)
		post = (alphaMat[1:,:] + betaMat[:-1,:]).transpose()
		post -= logsumexp(post,axis=0)
	return post

def likelihood(theta, args):
	Xvals = args[0]
	Yvals = args[1]
	numdimensions = round((np.sqrt( len(theta) * 8 + 1)  - 1)/2.0)
	if numdimensions == 1: #1d optimization
		standarddev = theta[0]
		scalarr  = 1/norm.pdf( args[2][0],  args[2][0], standarddev)
		if standarddev <=0:
			return(10000000000.0)
		FUNC = 0
		for i in range(len(Xvals)):
			FUNC = FUNC + (Yvals[i] - scalarr * norm.pdf(Xvals[i], loc = args[2][0], scale = standarddev))**2
		return(FUNC)
	else:		
		cholfactor = np.zeros((numdimensions, numdimensions))
		elementindex = 0

		for difference in range(numdimensions):
			for col in range(numdimensions - difference):
				cholfactor[col + difference, col] = theta[elementindex]
				elementindex = elementindex + 1
		for row in range(numdimensions):
			for col in range(numdimensions):
				if cholfactor[row, col] == 0.0:
					cholfactor[row, col] = cholfactor[col, row]
		if not np.all(np.linalg.eigvals(cholfactor) > 0):
			return(10000000000.0)
		try:
			scalarr  = 1/multivariate_normal.pdf(args[2], mean = args[2], cov = cholfactor)
		except:
			return(10000000000.0)

		FUNC = 0
		for i in range(len(Xvals)):
			try:
				FUNC = FUNC + (Yvals[i] - scalarr * multivariate_normal.pdf(Xvals[i], mean = args[2], cov = cholfactor))**2
			except:
				return(10000000000.0)
				
		return(FUNC)

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
	opts = {}

	if T == 2:
		Simplex = np.reshape(np.array([-0.05,0.05]),(2,1))
	elif T > 2:
		Simplex = np.zeros((T,T-1))
		for i in range(Simplex.shape[1]):
			Simplex[i,:] = 0.0
			Simplex[i,i] = 0.01
		Simplex[-1,:] = 0.0
	else:
		raise ValueError
	opts['disp']=False
	opts['maxfev'] = (T - 1) * 20
	opts['initial_simplex']=Simplex

	ImpSamp = False
	if times.shape[2] > 1:
		print('\t(Importance sampling with M = %d samples)'%(times.shape[2]))
		print()
		ImpSamp = True
	if not ImpSamp: # to account for whether we return the likelihood or the log likelihood

		logL0 = likelihood_wrapper(S0,timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals)

		minargs = (timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals)

		if len(S0) == 1:
			try:
				res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.9,1.0,1.1],args=minargs, method = "Brent", tol = 1e-4))
				S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
				L = res.fun

				if times.shape[1] < 31: # with small number of samples, result function can be multimodal, so do multiple optims and take best.
					try:
						res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [1.0001,1.00011,1.02],args=minargs, method = "Brent", tol = 1e-4))
						Stemp = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
						Ltemp = res.fun
						if Ltemp < L:
							S = Stemp
							L = Ltemp
					except ValueError:
						pass
			except ValueError:
				print("Selection MLE not found in [-0.1,0.1], possibly due to noninformative data. Expanding search to [-1,1].")
				res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.00000001,1.0,2.0],args=minargs, method = "Brent", tol = 1e-4))
				S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
				L = res.fun
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
		tranmatrix = _nstep_log_trans_prob(Ne[0],0.0,freqs,z_bins,z_logcdf,z_logsf,h) # this only handles the precompute case, just use initial values
		if len(np.unique(Ne)) == 1:
			precompute = 1
		for i in range(M):
			betaMatl0 = backward_algorithm(np.zeros(len(Ne)),times[:,:,i],derSampledTimes,ancSampledTimes,epochs,Ne,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,tranmatrix,noCoals=noCoals,precomputematrixboolean=precompute,currFreq=currFreq)
			Weights[i] = logsumexp(betaMatl0[-2,:])

		minargs = (timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,functionvals,Weights)

		if len(S0) == 1:
			try:
				res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.9,1.0,1.1],args=minargs, method = "Brent", tol = 1e-4))
				S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
				L = res.fun
				if times.shape[1] < 31: # with small number of samples, result function can be multimodal, so do multiple optims and take best.
					try:
						res = (minimize_scalar(likelihood_wrapper_scalar,  bracket = [1.0001,1.00011,1.02] ,args=minargs, method = "Brent", tol = 1e-4))
						Stemp = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
						Ltemp = res.fun
						if Ltemp < L:
							S = Stemp
							L = Ltemp
					except ValueError:
						pass
			except ValueError:
				print("Selection MLE not found in [-0.1,0.1], possibly due to noninformative data. Expanding search to [-1,1].")
				res = (minimize_scalar(likelihood_wrapper_scalar, bracket = [0.00000001,1.0,2.0],args=minargs, method = "Brent", tol = 1e-4))
				S = [res.x - 1.0] # adjusted wrapper to work on selection + 1, so that tolerance makes more sense.
				L = res.fun
		else:
			res = minimize(likelihood_wrapper, S0, args=minargs, options=opts, method='Nelder-Mead')
			S = res.x
			L = res.fun
		numericloglik = -L
		toprint = '%.4f'%(-L)

	FirstLine = "logLR" + "\t" + "-log10(p-value)"
	epochnum = 1
	degreesoffreedom = len(S)

	toprint = toprint + "\t" + '%.2f'%(-(chi2.logsf(numericloglik + numericloglik, degreesoffreedom ) ) / np.log(10) )
	for s,t,u in zip(S,timeBins[:-1],timeBins[1:]):
		toprint = toprint + "\t" + '%d'%(t)
		if u <= args.tCutoff:
			toprint = toprint + "\t" + '%d'%(u)
		else:
			toprint = toprint + "\t" + '%d'%(args.tCutoff)
		toprint = toprint + "\t" + '%.5f'%(s)
		FirstLine = FirstLine + "\t" + "Epoch" + str(epochnum) + "_start" + "\t" +  "Epoch" + str(epochnum)   + "_end"  + "\t" +  "SelectionMLE" + str(epochnum)
		epochnum = epochnum + 1
	f = open(args.out+"_inference.txt", "w+")
	f.writelines(FirstLine  + "\n" + toprint + "\n")
	f.close()
	functionvals.close()

	if not args.noAlleleTraj:
		file1 = open(args.out+"_tempfile.txt", 'r')
		Lines = file1.readlines()
		Xvals = []
		Yvals = []
		indexxx = 0
		for i in Lines:
			Xvals.append([])
			DerivedSampleTimes = i.split(",")
			for j in range(len(DerivedSampleTimes) - 1 ):
				(Xvals[indexxx]).append(float(DerivedSampleTimes[j]))
			Yvals.append(float(DerivedSampleTimes[len(DerivedSampleTimes) - 1 ]))
			indexxx = indexxx + 1
		Yvals = np.exp(np.subtract(Yvals, max(Yvals)))
		muu = Xvals[(list(Yvals)).index(max(Yvals))]

		if len(Xvals[0]) == 1:
			S0 =[1.0]
			res = minimize(likelihood, S0, args=[Xvals, Yvals, muu], method='Nelder-Mead', options={"maxfev":1000, "fatol":1e-20, "xatol":1e-20}).x
			#print("mu1: ", muu)
			#print("sd1: ", res[0])
			standard_dev = res[0]
			if args.integration_points == -1:
				variatessold = normal(loc=muu, scale=standard_dev, size=10)
			else:
				variatessold = normal(loc=muu, scale=standard_dev, size=args.integration_points)
			variatess = []
			for iiiv in variatessold:
				variatess.append([iiiv])
		else:
			S0 =[0.0] * (round((len(Xvals[0])*len(Xvals[0])+len(Xvals[0]))/2)  )
			for innn in range(len(Xvals[0]) ):
				S0[innn] = 1.0
			res = minimize(likelihood, S0, args=[Xvals, Yvals, muu], method='Nelder-Mead', options={"maxfev":1000, "fatol":1e-40, "xatol":1e-40})
			for ifi in range(1,9):
				S0 =[0.0] * (round((len(Xvals[0])*len(Xvals[0])+len(Xvals[0]))/2)  )
				for innn in range(len(Xvals[0])  ):
					S0[innn] = 10**(-ifi)
				res1 = minimize(likelihood, S0, args=[Xvals, Yvals, muu], method='Nelder-Mead', options={"maxfev":1000, "fatol":1e-40, "xatol":1e-40})

				if res1.fun < res.fun:
					res = res1
			#print("least-squares residual:", res.fun)
			res = res.x
			for iggi in res:
				if iggi > 0.2 or iggi < -0.2:
					print("Poor fit of normal distribution to data. Unreliable results follow.")

			#print("mu2: ", muu)
			#print("sd2: ", res)
			standard_dev = res
			numdimensions=  len(muu)

			covarmat = np.zeros((numdimensions, numdimensions))
			elementindex = 0

			for difference in range(numdimensions):
				for col in range(numdimensions - difference):
					covarmat[col + difference, col] = res[elementindex]
					elementindex = elementindex + 1
			for row in range(numdimensions):
				for col in range(numdimensions):
					if covarmat[row, col] == 0.0:
						covarmat[row, col] = covarmat[col, row]
			
			if args.integration_points == -1:
				variatess = multivariate_normal.rvs(mean=muu, cov=covarmat, size = numdimensions * 10)
			else:
				variatess = multivariate_normal.rvs(mean=muu, cov=covarmat, size = args.integration_points)

		# infer trajectory @ MLE of selection parameter
		post = np.exp(traj_wrapper(variatess[0],timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,Weights))
		for v in range(1,len(variatess)):
			post = post + np.exp(traj_wrapper(variatess[v],timeBins,Ne,h,freqs,times,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,sMax,derSampledTimes,ancSampledTimes,Weights))
		post = post / np.sum(post,axis=0)
		np.savetxt(args.out+"_freqs.txt", freqs, delimiter=",")
		np.savetxt(args.out+"_post.txt", post, delimiter=",")
	if os.path.exists(args.out+"_tempfile.txt"):
		os.remove(args.out+"_tempfile.txt")