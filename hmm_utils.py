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
	for i in range(lf):  ###can still significantly speed this up by maybe using a better truncator than 0.05 in either direction
		p1[i,:] = _log_trans_prob(i,N,s,FREQS,z_bins,z_logcdf,z_logsf)
	return(p1)

@njit('float64(float64[:],int64,float64[:],float64,float64,int64)',cache=True)
def _log_coal_density(times,n,epoch,xi,Ni,anc=0):
    if n == 1:
        # this flag indicates to ignore coalescence
        return 0.0
    #xi is current frequency, #Ni is DIPLOID pop size (1/2 of the previous one.)
    logp = 0
    prevt = epoch[0]
    if anc == 1:
        xi = 1.0-xi
    k=n
    xiNi = xi*Ni
    for i,t in enumerate(times):
        k = n-i
        kchoose2 = k*(k-1)/4
        dLambda = 1/xiNi*(t-prevt)
        logpk = - np.log(xi) - kchoose2 * dLambda
        logp += logpk
        
        prevt = t
        k -= 1
    kchoose2 = k*(k-1)/4
    logPk = - kchoose2 * 1/xiNi*(epoch[1]-prevt) # this is the survival function of the exponential. survival is exp(-log*x). chance of no further coalescences.

    logp += logPk
    return logp


@njit('float64[:,:](float64[:],float64[:,:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:],float64[:,:],float64[:,:],int64,int64,float64)',cache=True)
def backward_algorithm(sel,times,derSampledTimes,ancSampledTimes,epochs,N,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,TRANMATRIX, noCoals=1,precomputematrixboolean=0,currFreq=-1):

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
    alphaMat = np.full((T+1,lf), -1e20)
    alphaMat[0,:] = alpha
    
    prevNt = -1
    prevst = -1
    
    cumGens = 0
    
    nDerRemaining = np.sum(times[0,:]>=0) - len(derSampledTimes)
    nAncRemaining = np.sum(times[1,:]>=0) + 1 - len(ancSampledTimes)

    coalEmissions = np.zeros(lf)
    for tb in range(0,T):
        Nt = N[tb]
        epoch = np.array([cumGens,cumGens+1.0])
        st = sel[tb]
        prevAlpha = np.copy(alpha)
        if prevNt != Nt or prevst != st:
            if precomputematrixboolean:
                currTrans = TRANMATRIX
            else:
                currTrans = _nstep_log_trans_prob(Nt,st,freqs,z_bins,z_logcdf,z_logsf)   # note that the indexing here is reversed as opposed to the forward, because the usage of currtrans is transposed!!!
        
        #grab coal times during epoch
        # calculate coal emission probs
        if noCoals:
            coalEmissions = np.zeros(lf)
        else:
            derCoals = np.copy(times[0,:])
            derCoals = derCoals[derCoals > cumGens]
            derCoals = derCoals[derCoals <= cumGens+1.0]
            ancCoals = np.copy(times[1,:])
            ancCoals = ancCoals[ancCoals > cumGens]
            ancCoals = ancCoals[ancCoals <= cumGens+1.0]
            for j in range(lf):
                    if nDerRemaining > 1: # if multiple derived lineages, all is good
                        coalEmissions[j] = _log_coal_density(derCoals,nDerRemaining,epoch,freqs[j],Nt,anc=0)
                        coalEmissions[j] += _log_coal_density(ancCoals,nAncRemaining,epoch,freqs[j],Nt,anc=1)
                    elif nDerRemaining == 0 and nAncRemaining == 1:
                        if j != 0:
                            coalEmissions[j] = -1e20 # derived allele freq must be 0.
                    elif nDerRemaining == 0 and nAncRemaining > 1:
                        if j != 0:
                            coalEmissions[j] = -1e20 # derived allele freq must be 0.
                        else:
                            coalEmissions[j] = _log_coal_density(ancCoals,nAncRemaining,epoch,freqs[j],Nt,anc=1) # run only ancestral
                    elif nDerRemaining == 1 and nAncRemaining == 1:
                        if j != 0:
                            coalEmissions[j] = 0.0
                        else:
                            coalEmissions[j] = _log_coal_density(ancCoals,2,epoch,freqs[j],Nt,anc=1) # run with 2 ancestral lineages
                    elif nDerRemaining == 1 and nAncRemaining > 1:
                        if j != 0:
                            coalEmissions[j] = _log_coal_density(ancCoals,nAncRemaining,epoch,freqs[j],Nt,anc=1)
                        else:
                            coalEmissions[j] = _log_coal_density(ancCoals,nAncRemaining+1,epoch,freqs[j],Nt,anc=1) # run with 2 ancestral lineages
                    else:
                        print("Incorrect Polarization of Alleles!")
            
            numberofsampledder = derSampledTimes
            numberofsampledder = numberofsampledder[numberofsampledder > cumGens]
            numberofsampledder = len(numberofsampledder[numberofsampledder <= cumGens + 1.0])

            numberofsampledanc = ancSampledTimes
            numberofsampledanc = numberofsampledanc[numberofsampledanc > cumGens]
            numberofsampledanc = len(numberofsampledanc[numberofsampledanc <= cumGens + 1.0])

            #print(cumGens,numberofsampledder,numberofsampledanc)
            nDerRemaining -= len(derCoals)
            nAncRemaining -= len(ancCoals)

            nDerRemaining += numberofsampledder
            nAncRemaining += numberofsampledanc

        for i in range(0,lf):
            alpha[i] = _logsumexp(prevAlpha[0:lf] + currTrans[0:lf,i]) + coalEmissions[i]
            if np.isnan(alpha[i]):
                alpha[i] = -np.inf

        prevNt = Nt
        prevst = st
        cumGens += 1.0
        alphaMat[tb,:] = alpha
    return alphaMat