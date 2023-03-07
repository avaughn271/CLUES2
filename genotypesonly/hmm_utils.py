import numpy as np
from numba import njit

@njit('float64(float64[:])',cache=True)
def _logsumexp(a):
    """Standard logsumexp
    INPUT: a - vector of log probabilities
    OUTPUT: a number corresponding to log[exp(a1) + ... + exp(an)] but more stably computed
    """

    a_max = np.max(a)
    return np.log(np.sum(np.exp(a - a_max))) + a_max

@njit('float64(float64[:],float64[:])',cache=True)
def _logsumexpb(a,b):
    """Standard logsumexp
    INPUT: a - vector of log probabilities
           b - a vector of elements
    OUTPUT: a number corresponding to log[b1*exp(a1) + ... + bn*exp(an)] but more stably computed
    """
    a_max = np.max(a)
    return np.log(np.sum(b * np.exp(a - a_max))) + a_max

@njit('float64[:](int64,float64,float64,float64[:],float64[:],float64[:],float64[:])',cache=True)
def _log_trans_prob(i,N,s,FREQS,z_bins,z_logcdf,z_logsf):
    """Standard logsumexp
    INPUT: i - an index. ranges from 0 to df-1 inclusive. Index of frequency bin
           N - diploid effective popualtion size (smaller one)
           s - selection coefficient
           FREQS - frequency discrete values as computed earlier
           z_bins - some numbers spread between -1e5 and 1e5, more concentrated near 0. Length of 2198.
           z_logcdf - log of standard normal CDF of z_bins
           z_logsf - log of (1 - standard normal CDF) of z_bins
    OUTPUT: The set of len(FREQS) transition probabilities.
    conceptually, you calculate the resulting normal distribution. Any density greater than FREQS[df] is places on FREQS[df] 
    Any probability mass less than FREQS[0] is placed on FREQS[0] Otherwise, midpoints between frequency bins are used
    to round the resulting probability into a discrete bucket using the cdf.
    """
	# 1-generation transition prob based on Normal distn
	
    p = FREQS[i]    #starting frequency
    lf = len(FREQS)  # number of freq bins
    logP = np.NINF * np.ones(lf)  # initialize logP
    mu = p-s*p*(1.0-p) # this is the only place where dominance or selection changes come in 
    #            This should be the inverse of (x+sx)/(1+sx)
    # This is the mean of the normal distribution going back in time.
    sigma = np.sqrt(p*(1.0-p)/(4.0*N)) # variance of the normal approx

    logP[0] = np.interp(np.array([(FREQS[0]-mu)/sigma]),z_bins,z_logcdf)[0]
    logP[lf-1] = np.interp(np.array([(FREQS[lf-1]-mu)/sigma]),z_bins,z_logsf)[0]

    for j in range(1,lf-1):
        if j == 1:
            mlo = FREQS[0]
        else:
            mlo = np.mean(np.array([FREQS[j],FREQS[j-1]]))
        if j == lf-2:
            mhi = FREQS[j+1]
        else:
            mhi = np.mean(np.array([FREQS[j],FREQS[j+1]]))

        l1 = np.interp(np.array([(mlo-mu)/sigma]),z_bins,z_logcdf)[0]
        l2 = np.interp(np.array([(mhi-mu)/sigma]),z_bins,z_logcdf)[0]
        logP[j] = _logsumexpb(np.array([l1,l2]),np.array([-1.0,1.0]))

    return logP

@njit('float64[:,:](float64,float64,float64[:],float64[:],float64[:],float64[:])',cache=True)
def _nstep_log_trans_prob(N,s,FREQS,z_bins,z_logcdf,z_logsf):
	"""Same input as before, except i.
    Performs the above calculation on all possible input frequencies p.
    Output - A square matrix p1 of size df by df. Here p1[i,j] is the probability density
    backwards in time between the frequency index i to the frequency index j.
    """
	lf = len(FREQS)
	p1 = np.zeros((lf,lf))

	# load rows into p1
	for i in range(lf):
		p1[i,:] = _log_trans_prob(i,N,s,FREQS,z_bins,z_logcdf,z_logsf)

	return(p1)

@njit('float64(float64[:],float64)')
def _genotype_likelihood_emission(ancGLs,p):
	"""ancGLs is a list of size 3. p is the derived allele frequency.
	Returns the probability of the given emission"""
	logp = np.log(p)
	log1p = np.log(1-p)

	logGenoFreqs = np.array([log1p + log1p, 0.693147180559945309 + logp + log1p,logp + logp])
	emission = _logsumexp(logGenoFreqs + ancGLs)
	if np.isnan(emission):
		emission = -np.inf
	return emission

@njit('float64[:,:](float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:],float64[:,:])',cache=True)
def forward_algorithm(sel,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs):

    '''
    Moves forward in time from past to present
    '''
    
    lf = len(freqs)
    # neutral sfs
    alpha = np.ones(len(freqs))
    alpha -= _logsumexp(alpha)
    
    T = len(epochs)-1
    alphaMat = np.zeros((T,lf))
    alphaMat[-1,:] = alpha

    prevNt = -1
    prevst = -1
    
    cumGens = epochs[-1]
    for tb in range(T-2,-1,-1):
        Nt = N[tb]
        
        st = sel[tb]
        prevAlpha = np.copy(alpha)

        if prevNt != Nt or prevst != st:
            #change in selection/popsize, recalc trans prob
            currTrans = _nstep_log_trans_prob(Nt,st,freqs,z_bins,z_logcdf,z_logsf)
        
        #grab ancient GL rows
        ancientGLrows = ancientGLs[np.logical_and(ancientGLs[:,0] <= cumGens, ancientGLs[:,0] > cumGens - 1.0)]

        # calculate ancient GL emission probs
        glEmissions = np.zeros(lf)
        
        for j in range(lf):
            for iac in range(ancientGLrows.shape[0]):
                glEmissions[j] += _genotype_likelihood_emission(ancientGLrows[iac,1:],freqs[j])
                
        # calculate coal emission probs
        
        for i in range(lf):
            alpha[i] = _logsumexp(prevAlpha + currTrans[i,:] + glEmissions)
            if np.isnan(alpha[i]):
                alpha[i] = -np.inf
        
        prevNt = Nt
        prevst = st
        cumGens -= 1.0
        alphaMat[tb,:] = alpha
    return alphaMat
    
@njit('float64[:,:](float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:],float64)',cache=True)
def backward_algorithm(sel,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,currFreq=-1):

    '''
    Moves backward in time from present to past
    '''
    lf = len(freqs)
    alpha = np.zeros(lf)
    indexofcurrent = -1
    maxdistance = 2.0
    for i in range(lf):
        distt = abs(freqs[i] - currFreq) 
        alpha[i] = -1e20
        if distt < maxdistance:
             maxdistance = distt
             indexofcurrent = i
    alpha[indexofcurrent] = 0.0 # index of all -Infs except freq bin closest to the true current freq
    T = len(epochs)-1
    alphaMat = np.zeros((T,lf))
    
    prevNt = -1
    prevst = -1
    
    cumGens = 0
    
    for tb in range(0,T):
        Nt = N[tb]
        st = sel[tb]
        prevAlpha = np.copy(alpha)
        if prevNt != Nt or prevst != st:
            currTrans = _nstep_log_trans_prob(Nt,st,freqs,z_bins,z_logcdf,z_logsf)
        
        #grab ancient GL rows
        ancientGLrows = ancientGLs[ancientGLs[:,0] > cumGens]
        ancientGLrows = ancientGLrows[ancientGLrows[:,0] <= cumGens + 1.0]
        glEmissions = np.zeros(lf)
        for j in range(lf):
            for iac in range(ancientGLrows.shape[0]):
                glEmissions[j] += _genotype_likelihood_emission(ancientGLrows[iac,1:],freqs[j])

        for i in range(lf):
            alpha[i] = _logsumexp(prevAlpha + currTrans[:,i] ) + glEmissions[i]
            if np.isnan(alpha[i]):
                alpha[i] = -np.inf
        
        prevNt = Nt
        prevst = st
        
        cumGens += 1.0
        alphaMat[tb,:] = alpha
    return alphaMat