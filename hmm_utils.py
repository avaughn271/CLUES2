import numpy as np
from numba import njit

@njit('float64(float64[:])',cache=True)
def _logsumexp(a):
    a_max = np.max(a)
    return np.log(np.sum(np.exp(a - a_max))) + a_max

@njit('float64(float64[:],float64[:])',cache=True)
def _logsumexpb(a,b):
    a_max = np.max(a)
    return np.log(np.sum(b * np.exp(a - a_max))) + a_max

@njit('float64[:](int64,float64,float64,float64[:],float64[:],float64[:],float64[:])',cache=True)
def _log_trans_prob(i,N,s,FREQS,z_bins,z_logcdf,z_logsf):
	# 1-generation transition prob based on Normal distn
	
    p = FREQS[i]
    lf = len(FREQS)
    logP = np.NINF * np.ones(lf)
    mu = p-s*p*(1.0-p) # this is the only place where dominance or frequency changes come in 
    #            This should be the inverse of (x+sx)/(1+sx)
    sigma = np.sqrt(p*(1.0-p)/(4.0*N))

    pi0 = np.interp(np.array([(FREQS[0]-mu)/sigma]),z_bins,z_logcdf)[0]
    pi1 = np.interp(np.array([(FREQS[lf-1]-mu)/sigma]),z_bins,z_logsf)[0]

    middleP = np.zeros(lf-2)
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
        middleP[j-1] = _logsumexpb(np.array([l1,l2]),np.array([-1.0,1.0]))                    

    logP[0] = pi0
    logP[1:lf-1] = middleP
    logP[lf-1] = pi1

    return logP

@njit('float64[:,:](float64,float64,float64[:],float64[:],float64[:],float64[:])',cache=True)
def _nstep_log_trans_prob(N,s,FREQS,z_bins,z_logcdf,z_logsf):
	lf = len(FREQS)
	p1 = np.zeros((lf,lf))

	# load rows into p1
	for i in range(lf):
		p1[i,:] = _log_trans_prob(i,N,s,FREQS,z_bins,z_logcdf,z_logsf)

	return(p1)

@njit('float64(float64[:],float64)')
def _hap_genotype_likelihood_emission(ancGLs,p):
    logGenoFreqs = np.array([np.log(1-p),np.log(p)])
    emission = _logsumexp(logGenoFreqs + ancGLs)
    if np.isnan(emission):
        emission = -np.inf
    return emission

@njit('float64(float64[:],float64)')
def _genotype_likelihood_emission(ancGLs,p):
	logGenoFreqs = np.array([2*np.log(1-p),np.log(2) + np.log(p) + np.log(1-p),2*np.log(p)])
	emission = _logsumexp(logGenoFreqs + ancGLs)
	if np.isnan(emission):
		emission = -np.inf
	return emission

@njit('float64(float64[:],int64,float64[:],float64,float64,int64)',cache=True)
def _log_coal_density(times,n,epoch,xi,Ni,anc=0):
    if n == 1:
        # this flag indicates to ignore coalescence
        return 0.0
    
    logp = 0
    prevt = epoch[0]
    if anc == 1:
        xi = 1.0-xi
    k=n
    for i,t in enumerate(times):
        k = n-i
        kchoose2 = k*(k-1)/4
        dLambda = 1/(xi*Ni)*(t-prevt)
        logpk = - np.log(xi) - kchoose2 * dLambda
        logp += logpk
        
        prevt = t
        k -= 1
    kchoose2 = k*(k-1)/4
    logPk = - kchoose2 * 1/(xi*Ni)*(epoch[1]-prevt)

    logp += logPk
    return logp

@njit('float64[:,:](float64[:],float64[:,:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:],float64[:,:],int64)',cache=True)
def forward_algorithm(sel,times,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,noCoals=1):

    '''
    Moves forward in time from past to present
    '''
    
    lf = len(freqs)    
    # neutral sfs
    alpha = -np.log(freqs) 
    binEdges = np.array([0]+[0.5*(freqs[i]+freqs[i+1]) for i in range(len(freqs)-1)]+[1])
    alpha += np.log(np.diff(binEdges))
	
    alpha -= _logsumexp(alpha)
    
    T = len(epochs)-1
    alphaMat = np.zeros((T+1,lf))
    alphaMat[-1,:] = alpha

    prevNt = -1
    prevst = -1
    
    cumGens = epochs[-1]
    
    nDer = np.sum(times[0,:]>=0)+1
    nDerRemaining = nDer - np.sum(np.logical_and(times[0,:]>=0, times[0,:]<=epochs[-1]))
    nAnc = np.sum(times[1,:]>=0)+1
    nAncRemaining = nAnc - np.sum(np.logical_and(times[1,:]>=0, times[1,:]<=epochs[-1]))
    coalEmissions = np.zeros(lf)

    for tb in range(T-1,0,-1):
        epoch = np.array([cumGens - 1.0,cumGens])
        Nt = N[tb]
        
        st = sel[tb]
        prevAlpha = np.copy(alpha)

        if prevNt != Nt or prevst != st:
            #change in selection/popsize, recalc trans prob
            currTrans = _nstep_log_trans_prob(Nt,st,freqs,z_bins,z_logcdf,z_logsf)
        
        #grab ancient GL rows
        ancientGLrows = ancientGLs[np.logical_and(ancientGLs[:,0] <= cumGens, ancientGLs[:,0] > cumGens - 1.0)]
        ancientHapGLrows = ancientHapGLs[np.logical_and(ancientHapGLs[:,0] <= cumGens, ancientHapGLs[:,0] > cumGens - 1.0)]

        # calculate ancient GL emission probs
        glEmissions = np.zeros(lf)
        
        for j in range(lf):
            for iac in range(ancientGLrows.shape[0]):
                glEmissions[j] += _genotype_likelihood_emission(ancientGLrows[iac,1:],freqs[j])
            for iac in range(ancientHapGLrows.shape[0]):
                glEmissions[j] += _hap_genotype_likelihood_emission(ancientHapGLrows[iac,1:],freqs[j])
                
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
            nDerRemaining += len(derCoals)
            nAncRemaining += len(ancCoals)
            for j in range(lf):
                    coalEmissions[j] = _log_coal_density(derCoals,nDerRemaining,epoch,freqs[j],Nt,anc=0)
                    coalEmissions[j] += _log_coal_density(ancCoals,nAncRemaining,epoch,freqs[j],Nt,anc=1)

        for i in range(lf):
            alpha[i] = _logsumexp(prevAlpha + currTrans[i,:] + glEmissions + coalEmissions) 
            if np.isnan(alpha[i]):
                alpha[i] = -np.inf
        
        prevNt = Nt
        prevst = st
        cumGens -= 1.0
        alphaMat[tb,:] = alpha
    return alphaMat
    
@njit('float64[:,:](float64[:],float64[:,:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:],float64[:,:],int64,float64)',cache=True)
def backward_algorithm(sel,times,epochs,N,freqs,z_bins,z_logcdf,z_logsf,ancientGLs,ancientHapGLs,noCoals=1,currFreq=-1):

    '''
    Moves backward in time from present to past
    '''
    
    lf = len(freqs)
    alpha = np.zeros(lf)
    if currFreq != -1:
        nsamp = 1000
        for i in range(lf):
            k = int(currFreq*nsamp)
            alpha[i] = -np.sum(np.log(np.arange(2,k+1)))-np.sum(np.log(np.arange(2,nsamp-k+1)))+np.sum(np.log(np.arange(2,nsamp+1)))
            alpha[i] += k*np.log(freqs[i]) + (nsamp-k)*np.log(1-freqs[i])
            
    T = len(epochs)-1
    alphaMat = np.zeros((T+1,lf))
    alphaMat[0,:] = alpha
    
    prevNt = -1
    prevst = -1
    
    cumGens = 0
    
    nDer = np.sum(times[0,:]>=0)+1
    nDerRemaining = nDer
    nAnc = np.sum(times[1,:]>=0)+1
    nAncRemaining = nAnc
    coalEmissions = np.zeros(lf)

    for tb in range(0,T):
        Nt = N[tb]
        epoch = np.array([cumGens,cumGens+1.0])
        st = sel[tb]
        prevAlpha = np.copy(alpha)
        if prevNt != Nt or prevst != st:
            currTrans = _nstep_log_trans_prob(Nt,st,freqs,z_bins,z_logcdf,z_logsf)
        
        #grab ancient GL rows
        ancientGLrows = ancientGLs[ancientGLs[:,0] > cumGens]
        ancientGLrows = ancientGLrows[ancientGLrows[:,0] <= cumGens + 1.0]

        ancientHapGLrows = ancientHapGLs[ancientHapGLs[:,0] > cumGens]
        ancientHapGLrows = ancientHapGLrows[ancientHapGLrows[:,0] <= cumGens + 1.0]

        glEmissions = np.zeros(lf)
        for j in range(lf):
            for iac in range(ancientGLrows.shape[0]):
                glEmissions[j] += _genotype_likelihood_emission(ancientGLrows[iac,1:],freqs[j])
            for iac in range(ancientHapGLrows.shape[0]):
                glEmissions[j] += _hap_genotype_likelihood_emission(ancientHapGLrows[iac,1:],freqs[j])
        
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
                    coalEmissions[j] = _log_coal_density(derCoals,nDerRemaining,epoch,freqs[j],Nt,anc=0)
                    coalEmissions[j] += _log_coal_density(ancCoals,nAncRemaining,epoch,freqs[j],Nt,anc=1)
            nDerRemaining -= len(derCoals)
            nAncRemaining -= len(ancCoals)

        for i in range(lf):
            alpha[i] = _logsumexp(prevAlpha + currTrans[:,i] ) + glEmissions[i] + coalEmissions[i]
            if np.isnan(alpha[i]):
                alpha[i] = -np.inf
        
        prevNt = Nt
        prevst = st
        
        cumGens += 1.0
        alphaMat[tb,:] = alpha
    return alphaMat

@njit('float64(float64[:,:],float64[:],float64[:])',cache=True)
def proposal_density(times,epochs,N):
    '''
    Moves backward in time from present to past
    '''
    
    logl = 0.
    cumGens = 0
    combinedTimes = np.sort(np.concatenate((times[0,:],times[1,:])))
    nRemaining = np.sum(combinedTimes>=0) + 2
    for tb in range(0,len(epochs)-1):
        Nt = N[tb]
        epoch = np.array([cumGens,cumGens+1.0])
        
        #grab coal times during epoch
        # calculate coal emission probs

        Coals = np.copy(combinedTimes)
        Coals = Coals[Coals > cumGens]
        Coals = Coals[Coals <= cumGens+1.0]
      
        logl += _log_coal_density(Coals,nRemaining,epoch,1.0,Nt,anc=0)
        nRemaining -= len(Coals)
        
        cumGens += 1.0
    return logl
