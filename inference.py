import numpy as np
from hmm_utils import forward_algorithm
from hmm_utils import backward_algorithm
from hmm_utils import proposal_density
from scipy.special import logsumexp
from scipy.optimize import minimize
import argparse
import gzip
import os

def parse_clues(filename,args):
    with gzip.open(filename, 'rb') as fp:
        try:
            #parse file
            data = fp.read()
        except OSError:
            with open(filename, 'rb') as fp:
                try:
                    #parse file
                    data = fp.read()
                except OSError:
                    print('Error: Unable to open ' + filename)
                    exit(1)
           
        #get #mutations and #sampled trees per mutation
        filepos = 0
        num_muts, num_sampled_trees_per_mut = np.frombuffer(data[slice(filepos, filepos+8, 1)], dtype = np.int32)

        filepos += 8
        #iterate over mutations
        for m in range(0,num_muts):
            bp = np.frombuffer(data[slice(filepos, filepos+4, 1)], dtype = np.int32)
            filepos += 4
            anc, der = np.frombuffer(data[slice(filepos, filepos+2, 1)], dtype = 'c')
            filepos += 2
            daf, n = np.frombuffer(data[slice(filepos, filepos+8, 1)], dtype = np.int32)
            filepos += 8
            
            if daf >= n-1:
            	anctimes = np.empty((num_sampled_trees_per_mut,0))
            else:
                num_anctimes = 4*(n-daf-1)*num_sampled_trees_per_mut
                anctimes     = np.reshape(np.frombuffer(data[slice(filepos, filepos+num_anctimes, 1)], dtype = np.float32), (num_sampled_trees_per_mut, n-daf-1))
                filepos     += num_anctimes
            
            if daf <= 1:
            	dertimes = np.empty((num_sampled_trees_per_mut,0))
            else:
                num_dertimes = 4*(daf-1)*num_sampled_trees_per_mut
                dertimes     = np.reshape(np.frombuffer(data[slice(filepos, filepos+num_dertimes, 1)], dtype = np.float32), (num_sampled_trees_per_mut, daf-1))
                filepos     += num_dertimes
		
    return dertimes,anctimes

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--times',type=str,
		help='Should refer to files <times>.{{der,anc}}.npy (exclude prefix .{{der,anc}}.npy)',
		default=None)
	parser.add_argument('--popFreq',type=float,default=None)

	parser.add_argument('--ancientSamps',type=str,default=None)
	parser.add_argument('--ancientHaps',type=str,default=None)
	parser.add_argument('--out',type=str,default=None)

	parser.add_argument('-N','--N',type=float,default=10**4)
	parser.add_argument('-coal','--coal',type=str,default=None,help='path to Relate .coal file. Negates --N option.')
	parser.add_argument('--dom',type=float,default=0.5,help='dominance coefficient')

	# adv options
	parser.add_argument('-thin','--thin',type=int,default=1)
	parser.add_argument('-burnin','--burnin',type=int,default=0)
	parser.add_argument('--tCutoff',type=float,default=1000)
	parser.add_argument('--timeBins',type=str,default=None)
	parser.add_argument('--sMax',type=float,default=0.1)
	parser.add_argument('--df',type=int,default=150)
	return parser.parse_args()
    

def load_times(args):
	locusDerTimes,locusAncTimes = parse_clues(args.times+'.timeb',args)
	print(locusDerTimes.shape,locusAncTimes.shape)	
	if locusDerTimes.ndim == 0 or locusAncTimes.ndim == 0:
		raise ValueError
	elif locusAncTimes.ndim == 1 and locusDerTimes.ndim == 1:
		M = 1
		locusDerTimes = np.transpose(np.array([locusDerTimes]))
		locusAncTimes = np.transpose(np.array([locusAncTimes]))
	elif locusAncTimes.ndim == 2 and locusDerTimes.ndim == 1:
		locusDerTimes = np.array([locusDerTimes])[:,args.burnin::args.thin]
		locusAncTimes = np.transpose(locusAncTimes)[:,args.burnin::args.thin]
		M = locusDerTimes.shape[1]	
	elif locusAncTimes.ndim == 1 and locusDerTimes.ndim == 2:
		locusAncTimes = np.array([locusAncTimes])[:,args.burnin::args.thin]
		locusDerTimes = np.transpose(locusDerTimes)[:,args.burnin::args.thin]
		M = locusDerTimes.shape[1]
	else:
		locusDerTimes = np.transpose(locusDerTimes)[:,args.burnin::args.thin]
		locusAncTimes = np.transpose(locusAncTimes)[:,args.burnin::args.thin]
		M = locusDerTimes.shape[1]
	n = locusDerTimes.shape[0] + 1
	m = locusAncTimes.shape[0] + 1
	ntot = n + m
	row0 = -1.0 * np.ones((ntot,M))


	row0[:locusDerTimes.shape[0],:] = locusDerTimes

	row1 = -1.0 * np.ones((ntot,M))

	row1[:locusAncTimes.shape[0],:] = locusAncTimes
	locusTimes = np.array([row0,row1])
	return locusTimes

def load_data(args):
		# load coalescence times
	noCoals = (args.times == None)
	if not noCoals:
		times = load_times(args)
	else:
		times = np.zeros((2,0,0))
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

	if noCoals:
		try:
			tCutoff = np.max(ancientGLs[:,0])+1.0
		except:
			tCutoff = np.max(ancientHapGLs[:,0])+1.0
	else:
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
	c = 1/(2*np.min([Ne[0],100000]))
	df = args.df
	freqs = np.linspace(c,1-c,df)
	# load time bins (for defining selection epochs)
	if args.timeBins != None:
		timeBins = np.genfromtxt(args.timeBins)
	else:
		timeBins = np.array([0.0,tCutoff])

	return timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,args.dom

def likelihood_wrapper(theta,timeBins,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,h,sMax):
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
    	for i in range(M):
    		betaMat = backward_algorithm(sel,times[:,:,i],epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals,currFreq=currFreq,h=h)
    		logl = logsumexp(betaMat[-2,:])
    		logl0 = proposal_density(times[:,:,i],epochs,N)
    		loglrs[i] = logl-logl0
    	logl = -1 * (-np.log(M) + logsumexp(loglrs))
    else:
    	betaMat = backward_algorithm(sel,t,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals,currFreq=currFreq,h=h)
    	logl = -logsumexp(betaMat[-2,:])
    return logl

def traj_wrapper(theta,timeBins,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,gens,noCoals,currFreq,h,sMax):
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
    	for i in range(M):
    		betaMat = backward_algorithm(sel,times[:,:,i],epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,np.array([]),noCoals=noCoals,currFreq=currFreq,h=h)
    		alphaMat = forward_algorithm(sel,times[:,:,i],epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals,h=h)
    		logl = logsumexp(betaMat[-2,:])
    		logl0 = proposal_density(times[:,:,i],epochs,N)
    		loglrs[i] = logl-logl0
    		postBySamples[:,:,i] = (alphaMat[1:,:] + betaMat[:-1,:]).transpose()
    	post = logsumexp(loglrs + postBySamples,axis=2)
    	post -= logsumexp(post,axis=0)

    else:
    	post = np.zeros((F,T))
    	betaMat = backward_algorithm(sel,t,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals,currFreq=currFreq,h=h)
    	alphaMat = forward_algorithm(sel,t,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancGLs,ancHapGLs,noCoals=noCoals,h=h)
    	post = (alphaMat[1:,:] + betaMat[:-1,:]).transpose()
    	post -= logsumexp(post,axis=0)
    return post

if __name__ == "__main__":
	args = parse_args()
	if args.times == None and args.ancientSamps == None and args.ancientHaps == None:
		print('You need to supply coalescence times (--times) and/or ancient samples (--ancientSamps) and/or ancient haploid samples (--ancientHaps)')
	
	print()
	print('Loading data and initializing model...')

	# load data and set up model
	sMax = args.sMax	
	timeBins,times,epochs,Ne,freqs,ancientGLs,ancientHapGLs,noCoals,currFreq,h = load_data(args)
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
	    
	logL0 = likelihood_wrapper(S0,timeBins,Ne,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,h,sMax)

	print('Optimizing likelihood surface using Nelder-Mead...')
	if times.shape[2] > 1:
		print('\t(Importance sampling with M = %d Relate samples)'%(times.shape[2]))
		print()
	minargs = (timeBins,Ne,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,h,sMax)
	res = minimize(likelihood_wrapper, S0, args=minargs, options=opts, method='Nelder-Mead')

	S = res.x
	L = res.fun

	print('#'*10)
	print()
	print('logLR: %.4f'%(-res.fun+logL0))
	print()
	print('MLE:')
	print('========')
	print('epoch\tselection')
	for s,t,u in zip(S,timeBins[:-1],timeBins[1:]):
		print('%d-%d\t%.5f'%(t,u,s))

	# infer trajectory @ MLE of selection parameter
	print(noCoals)

	post = traj_wrapper(res.x,timeBins,Ne,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,epochs,noCoals,currFreq,h,sMax)
	
	if args.out != None:
		np.save(args.out+'.epochs',epochs)
		np.save(args.out+'.freqs',freqs)
		np.save(args.out+'.post',post)
	else:
		print()
		print('Trajectory:')
		print('=============')
		print('gens_bp\tfreq')
		for i in range(0,int(timeBins[-1]),int(timeBins[-1]//(50))):
			if i > epochs[-1]:
				break
			print(i,np.sum(freqs * np.exp(post[:,i])))
		print()
		print('Finished.')
		print()
