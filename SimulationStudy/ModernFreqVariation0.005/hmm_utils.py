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

@njit('float64(float64,float64,float64,float64[:],float64[:])',cache=True)
def general_normal_cdf(x,mean, sd, xvals, yvals):
    return np.interp( (x - mean)/sd , xvals, yvals)

@njit('float64[:](float64[:],int64,float64,float64,float64[:],float64[:],float64[:],float64[:],float64)',cache=True)
def _log_trans_prob(BINGAPS, i,N,s,FREQS,z_bins,z_logcdf,z_logsf,h):
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
    logP = 0.0 * np.ones(lf)  # initialize logP, not log at the beginning.
    #mu = p-s*p*(1.0-p) # this is the only place where dominance or selection changes come in 
    #            This should be the inverse of (x+sx)/(1+sx)
    #print("h is:", h)
    mu = p + (s* (-1 + p)* p* (-p + h *(-1 + 2* p))) /(-1 + s *(2* h *(-1 + p) - p) *p) # general dominance coefficient
    # This is the mean of the normal distribution going back in time.
    sigma = np.sqrt(p*(1.0-p)/(2.0*N)) #IS THIS RIGHT????? maybe change back to 4 * N

    lowerfrequencybound = mu - sigma * 3.3 # guarantees 99.9 of probability is computed.
    upperfrequnecybound = mu + sigma * 3.3
    lowerindex = np.argmin(np.abs(np.subtract(FREQS, lowerfrequencybound))) - 1
    upperindec = np.argmin(np.abs(np.subtract(FREQS, upperfrequnecybound))) + 1
    upperindec = min(upperindec, lf)
    lowerindex = max(lowerindex, 0)

    for j in range(lowerindex,upperindec):
        if j == 0:
            logP[j] = general_normal_cdf(BINGAPS[j],   mu,   sigma, z_bins, z_logcdf)
        elif j == lf - 1:
            logP[lf-1] = 1 - general_normal_cdf(BINGAPS[-1],   mu,   sigma, z_bins, z_logcdf)
        else:
            logP[j] = general_normal_cdf(BINGAPS[j], mu, sigma, z_bins, z_logcdf) - general_normal_cdf(BINGAPS[j - 1], mu, sigma, z_bins, z_logcdf)
    return np.log(logP / np.sum(logP)) # renormalize, change so we only compute log of relevant entries??

@njit('float64[:,:](float64,float64,float64[:],float64[:],float64[:],float64[:],float64)',cache=True)
def _nstep_log_trans_prob(N,s,FREQS,z_bins,z_logcdf,z_logsf,h):
	"""Same input as before, except i.
    Performs the above calculation on all possible input frequencies p.
    Output - A square matrix p1 of size df by df. Here p1[i,j] is the probability density
    backwards in time between the frequency index i to the frequency index j.
    """
	lf = len(FREQS)
	p1 = np.zeros((lf,lf))
    
	BINGAPS = (FREQS[1:] + FREQS[0:(len(FREQS) - 1)]) / 2.0
	p1[0,0] = 1.0 # these are obsorbing states.
	p1[lf-1,lf-1] = 1.0
	p1[0,:] = np.log(p1[0,:])
	p1[lf-1,:] = np.log(p1[lf-1,:])

	# load rows into p1
	for i in range(1, lf - 1):  ###can still significantly speed this up by maybe using a better truncator than 0.05 in either direction
		p1[i,:] = _log_trans_prob(BINGAPS, i,N,s,FREQS,z_bins,z_logcdf,z_logsf,h)

	return(p1)

@njit('float64(float64[:],float64,float64)')
def _hap_genotype_likelihood_emission(ancGLs,logp, log1p):
	"""Be aware that this does not exactly convert to a corresponding genotype emission due to specifying which haplotype is which. A factor of 2 is the difference."""
	logGenoFreqs = np.array([log1p,logp])
	emission = _logsumexp(logGenoFreqs + ancGLs)
	if np.isnan(emission):
		emission = -np.inf
	return emission

@njit('float64(float64[:],float64,float64)')
def _genotype_likelihood_emission(ancGLs,logp, log1p):
	"""ancGLs is a list of size 3. p is the derived allele frequency.
	Returns the probability of the given emission"""

	logGenoFreqs = np.array([log1p + log1p, 0.693147180559945309 + logp + log1p,logp + logp])
	emission = _logsumexp(logGenoFreqs + ancGLs)
	if np.isnan(emission):
		emission = -np.inf
	return emission

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
    if xiNi == 0.0: # if in other state, impossible for others to still be segregating
        return(-1e20)
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

@njit('float64[:,:](float64[:],float64[:,:],float64[:],float64[:],float64[:],float64[:],float64,float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:],float64[:,:],int64)',cache=True)
def forward_algorithm(sel,times,derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,noCoals=1):

    '''
    Moves forward in time from past to present
    ''' 
    
    lf = len(freqs)
    alpha = np.ones(len(freqs)) # uniform probability of exiting at furthest time point.
    alpha -= _logsumexp(alpha)
    
    T = len(epochs)-1
    alphaMat = np.full((T+1,lf), -1e20)
    alphaMat[-1,:] = alpha

    prevNt = -1
    prevst = -1
    
    cumGens = epochs[-1]
    
    nDer = np.sum(times[0,:]>=0)
    nDerRemaining = nDer - np.sum(np.logical_and(times[0,:]>=0, times[0,:]<=epochs[-1]))  - len(derSampledTimes[derSampledTimes > epochs[-1]])          ###??????????????
    nAnc = np.sum(times[1,:]>=0)+1
    nAncRemaining = nAnc - np.sum(np.logical_and(times[1,:]>=0, times[1,:]<=epochs[-1]))  - len(ancSampledTimes[ancSampledTimes > epochs[-1]])  #OR PLUS??   
    coalEmissions = np.zeros(lf)

    for tb in range(T-1,0,-1):
        epoch = np.array([cumGens - 1.0,cumGens])
        Nt = N[tb]
        
        st = sel[tb]
        prevAlpha = np.copy(alpha)

        if prevNt != Nt or prevst != st:
            #change in selection/popsize, recalc trans prob
            currTrans = _nstep_log_trans_prob(Nt,st,freqs,z_bins,z_logcdf,z_logsf,h)
            upperindex = np.zeros(lf)
            lowerindex = np.zeros(lf)
            neginf = np.log(0.0)
            for ii in range(lf):  #calculate the range of the nonzero entries of the transition matrix.
                chosenrow =  currTrans[ii,:]
                lowerbound = np.argmax(currTrans[ii,:])
                upperbound = lowerbound + 1
                while (chosenrow[lowerbound] > neginf):
                    if lowerbound == 0:
                        break
                    lowerbound = lowerbound - 1
                while (chosenrow[upperbound] > neginf):
                    if upperbound == lf:
                        break
                    upperbound = upperbound + 1

                upperindex[ii] = upperbound
                lowerindex[ii] = lowerbound
        
        #grab ancient GL rows
        ancientGLrows = ancientGLs[np.logical_and(ancientGLs[:,0] <= cumGens, ancientGLs[:,0] > cumGens - 1.0)]
        ancientHapGLrows = ancientHapGLs[np.logical_and(ancientHapGLs[:,0] <= cumGens, ancientHapGLs[:,0] > cumGens - 1.0)]


        ancientHapGLrowsDerived = derSampledTimes[derSampledTimes > cumGens - 1.0]
        ancientHapGLrowsDerived = ancientHapGLrowsDerived[ancientHapGLrowsDerived <= cumGens]
        ancientHapGLrowsAncestral = ancSampledTimes[ancSampledTimes > cumGens - 1.0]
        ancientHapGLrowsAncestral = ancientHapGLrowsAncestral[ancientHapGLrowsAncestral <= cumGens]



        # calculate ancient GL emission probs
        glEmissions = np.zeros(lf)
        
        for j in range(lf):
            for iac in range(ancientGLrows.shape[0]):
                glEmissions[j] += _genotype_likelihood_emission(ancientGLrows[iac,1:],logfreqs[j],log1minusfreqs[j])
            for iac in range(ancientHapGLrows.shape[0]):
                glEmissions[j] += _hap_genotype_likelihood_emission(ancientHapGLrows[iac,1:],logfreqs[j],log1minusfreqs[j])
            for iac in range(len(ancientHapGLrowsDerived)):
                glEmissions[j] += np.log(freqs[j])
            for iac in range(len(ancientHapGLrowsAncestral)):
                glEmissions[j] += np.log(1 - freqs[j])
                
        # calculate coal emission probs
        
        if noCoals:
            coalEmissions = np.zeros(lf)
        else:
            derCoals = np.copy(times[0,:])
            derCoals = derCoals[derCoals <= cumGens]
            derCoals = derCoals[derCoals > cumGens-1.0]
            ancCoals = np.copy(times[1,:])
            ancCoals = ancCoals[ancCoals <= cumGens]
            ancCoals = ancCoals[ancCoals > cumGens-1.0]



            numberofsampledder = derSampledTimes
            numberofsampledder = numberofsampledder[numberofsampledder > cumGens - 1.0]

            numberofsampledanc = ancSampledTimes
            numberofsampledanc = numberofsampledanc[numberofsampledanc > cumGens - 1.0]

            nDerRemaining += len(derCoals)
            nAncRemaining += len(ancCoals)

            nDerRemaining -= len(numberofsampledder[numberofsampledder <= cumGens])
            nAncRemaining -= len(numberofsampledanc[numberofsampledanc <= cumGens])

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
        for i in range(lf):
            alpha[i] = _logsumexp(prevAlpha[lowerindex[i]:upperindex[i]] + currTrans[i,lowerindex[i]:upperindex[i]] + glEmissions[lowerindex[i]:upperindex[i]] + coalEmissions[lowerindex[i]:upperindex[i]])
            if np.isnan(alpha[i]):
                alpha[i] = -np.inf
        
        prevNt = Nt
        prevst = st
        cumGens -= 1.0
        alphaMat[tb,:] = alpha
    return alphaMat

@njit('float64[:,:](float64[:],float64[:,:],float64[:],float64[:],float64[:],float64[:],float64,float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:],float64[:,:],float64[:,:],int64,int64,float64)',cache=True)
def backward_algorithm(sel,times,derSampledTimes,ancSampledTimes,epochs,N,h,freqs,logfreqs,log1minusfreqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,TRANMATRIX, noCoals=1,precomputematrixboolean=0,currFreq=-1):

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
                currTrans = _nstep_log_trans_prob(Nt,st,freqs,z_bins,z_logcdf,z_logsf,h)   # note that the indexing here is reversed as opposed to the forward, because the usage of currtrans is transposed!!!

            upperindex = np.zeros(lf)
            lowerindex = np.zeros(lf)
            neginf = np.log(0.0)
            for ii in range(lf):  #calculate the range of the nonzero entries of the transition matrix.
                chosenrow =  currTrans[:,ii]
                lowerbound = np.argmax(currTrans[:,ii])
                upperbound = lowerbound + 1
                while (chosenrow[lowerbound] > neginf):
                    if lowerbound == 0:
                        break
                    lowerbound = lowerbound - 1
                if upperbound != lf: # added
                    while (chosenrow[upperbound] > neginf):
                        if upperbound == lf:
                            break
                        upperbound = upperbound + 1
                upperindex[ii] = upperbound
                lowerindex[ii] = lowerbound
        #grab ancient GL rows
        ancientGLrows = ancientGLs[ancientGLs[:,0] > cumGens]
        ancientGLrows = ancientGLrows[ancientGLrows[:,0] <= cumGens + 1.0]
        ancientHapGLrows = ancientHapGLs[ancientHapGLs[:,0] > cumGens]
        ancientHapGLrows = ancientHapGLrows[ancientHapGLrows[:,0] <= cumGens + 1.0]
        ancientHapGLrowsDerived = derSampledTimes[derSampledTimes > cumGens]
        ancientHapGLrowsDerived = ancientHapGLrowsDerived[ancientHapGLrowsDerived <= cumGens + 1.0]
        ancientHapGLrowsAncestral = ancSampledTimes[ancSampledTimes > cumGens]
        ancientHapGLrowsAncestral = ancientHapGLrowsAncestral[ancientHapGLrowsAncestral <= cumGens + 1.0]
        glEmissions = np.zeros(lf)
        for j in range(lf):
            for iac in range(ancientGLrows.shape[0]):
                glEmissions[j] += _genotype_likelihood_emission(ancientGLrows[iac,1:],logfreqs[j],log1minusfreqs[j])
            for iac in range(ancientHapGLrows.shape[0]):
                glEmissions[j] += _hap_genotype_likelihood_emission(ancientHapGLrows[iac,1:],logfreqs[j],log1minusfreqs[j])
            for iac in range(len(ancientHapGLrowsDerived)):
                glEmissions[j] += np.log(freqs[j])
            for iac in range(len(ancientHapGLrowsAncestral)):
                glEmissions[j] += np.log(1 - freqs[j])
        
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

            numberofsampledanc = ancSampledTimes
            numberofsampledanc = numberofsampledanc[numberofsampledanc > cumGens]

            nDerRemaining -= len(derCoals)
            nAncRemaining -= len(ancCoals)

            nDerRemaining += len(numberofsampledder[numberofsampledder <= cumGens + 1.0])
            nAncRemaining += len(numberofsampledanc[numberofsampledanc <= cumGens + 1.0])
        if h > 0.4 and h < 0.6:
            maxloc = np.argmax(prevAlpha)
            previouscolumn = np.exp(prevAlpha - prevAlpha[maxloc])
            prevsumtarget = np.sum(previouscolumn) * 0.999
            runningsum = previouscolumn[maxloc]
            currlower = maxloc - 1
            currhigher = maxloc + 1
            while (runningsum < prevsumtarget):
                if currhigher > lf - 1:
                    runningsum = runningsum + previouscolumn[currlower]
                    currlower = currlower - 1
                elif currlower < 0:
                    runningsum = runningsum + previouscolumn[currhigher]
                    currhigher = currhigher + 1
                else:
                    lowerpossible =  previouscolumn[currlower]
                    upperpossible =  previouscolumn[currhigher]
                    if lowerpossible < upperpossible:
                        runningsum = runningsum + previouscolumn[currhigher]
                        currhigher = currhigher + 1
                    else:
                        runningsum = runningsum + previouscolumn[currlower]
                        currlower = currlower - 1
            if currlower < 0:
                currlower = currlower + 1
            if currhigher > lf -1:
                currhigher = currhigher - 1
            lowerbounddd = 0
            upperbounddd = lf
            if freqs[currhigher] - freqs[currlower] < 0.33:
                indexofcurrent = -1
                maxdistance = 100.0
                target =  freqs[currlower]  - 0.1 - 2*freqs[currhigher]  + 2*freqs[currlower]
                for i in range(lf):
                    distt = abs(freqs[i] - target) 
                    if distt < maxdistance:
                        maxdistance = distt
                        indexofcurrent = i
                lowerbounddd = indexofcurrent

                indexofcurrent = -1
                maxdistance = 100.0
                target =  freqs[currhigher] + 0.1 + 2*freqs[currhigher] -  2*freqs[currlower]
                for i in range(lf):
                    distt = abs(freqs[i] - target)
                    if distt < maxdistance:
                        maxdistance = distt
                        indexofcurrent = i
                upperbounddd = indexofcurrent
                upperbounddd = upperbounddd + 1

            if nDerRemaining == 1:
                lowerbounddd = min(lowerbounddd, 0)
        else:
            lowerbounddd = 0
            upperbounddd = lf
        ####################################
        for i in range(lowerbounddd,upperbounddd):
            alpha[i] = _logsumexp(prevAlpha[lowerindex[i]:upperindex[i]] + currTrans[lowerindex[i]:upperindex[i],i]) + glEmissions[i] + coalEmissions[i]
            if np.isnan(alpha[i]):
                alpha[i] = -np.inf

        prevNt = Nt
        prevst = st
        cumGens += 1.0
        alphaMat[tb,:] = alpha
    return alphaMat
