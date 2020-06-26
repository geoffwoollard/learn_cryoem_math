#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Jun  12 2019
@author: pcn  cryo1dMaxLik.py  Follow Sigworth J Str Biol 1998 but in 1D
Python 3.7
Description: This code does not attempt to refine estimates of drift nor noise strength; they are taken as given.
create image 1dbests
"""
import numpy as np; import matplotlib.pyplot as plt
from scipy import ndimage

# parameters
iexamine = 1 # choose one image to show details about
Nrefine = 16

plt.close('all')
a = np.load('1dImages.npz')
samples = a['samples']   # fake datasets [which image in dataset, which pixel in image, which dataset]
Nsamp,siz,Nnois = samples.shape
if siz%2==0: raise ValueError('siz should be odd')

shiftSD = a['shiftSD'] # SD of the jitter in each fake dataset
noiselevels = a['noiselevels'] # 1/sqrt(SNR) in each pixel of each dataset
noiseSDs = noiselevels/np.sqrt(1+noiselevels**2)


def reestimate(prevTheta, images):
    """refines prevTheta=[A,noisesig,xisig,xiq] and returns newTheta and loglik value
    images[whicch image, x] is a set of data samples for a particular noise level, e.g. samples[:,:,1]
    this simplified code does not attempt to update noisesig, xisig nor xiq (Sigworth paper eq 17-19)
"""
    Aprev = prevTheta[0] # previous estimated image
    siz = len(Aprev); halfsiz = int((siz-1)/2)
    Nsamp, sizb = images.shape
    if sizb!=siz: raise ValueError('size mismatch')
    
    noisesig = prevTheta[1] # previous estimate of noise 
    xisig = prevTheta[2] # previous estimate of position jitter
    xiq = prevTheta[3] # previous estimate of position centering
    
    xisigalt = 1./(2.*xisig**2) # convenient
    noisesigalt = 1./noisesig**2
    imagenormtot = np.sum(images*images) # sum of norms of all images
    
    # tabulate Blook =  norms of shifted targets
    Blookup = np.zeros(siz)
    Asquares = Aprev*Aprev
    Blookup[halfsiz] = np.sum(Asquares) # the biggest norm is the unshifted one
    for sh in range(halfsiz): # could vectorize with cumsum
        Blookup[sh] =  np.sum(Asquares[:(halfsiz+sh+1)]) #  shifts (rightward)
        Blookup[siz-1-sh] = np.sum(Asquares[(halfsiz-sh):]) # shifts (leftward) 
    Blook = Blookup.reshape((1,siz)) # convert to a row vector
    xs = np.arange(siz)
    xss = np.reshape(xs, (1,siz)) # convert to a row vector
    gambarB = np.exp(-xisigalt*(xss-halfsiz-xiq)**2) # broadcast (ind of which image i)
    
    # correlate every micrograph (=images, the simulated exp data) with Aprev:
    corsX = ndimage.correlate1d(images, Aprev, axis=1, mode='constant', cval=0.)
    # actually we want Corr(Aprev, images) but correlate1d won't allow that directly so fix:
    cors = corsX[:,-1::-1]
#    Ki = noisesigalt*(-Blookup[halfsiz]/2 + cors[:,halfsiz]) #center values, broadcast 1st term
    Ki = np.amax(noisesigalt*(-Blook/2 + cors), axis=1) # refer to this common value
    gambar = np.exp(noisesigalt*(-Blook/2 + cors) - Ki.reshape((Nsamp,1)))*gambarB # broadcast 1st and 3rd terms in exp
    gambarsum = np.sum(gambar, axis=1)
    # compute new loglik: reinstate the terms that were subtracted/omitted:
    LL = np.sum(np.log(gambarsum) + Ki) - noisesigalt*imagenormtot/2
    # DEBUG:
    print ('best gammabar for image ',iexamine,' is at ',np.argmax(gambar[iexamine,:])) # show best shift
    # compute updated image:
    Anext = np.zeros(siz)
    for i in range(Nsamp):
        Anext += ndimage.convolve1d(gambar[i,:], images[i,:])/gambarsum[i]
    Anext = Anext/Nsamp
    return ([Anext, noisesig, xisig, xiq], LL)  # temporarily don't update last 3

"""main body (calls reestimate)"""
avgs = np.mean(samples, axis=0) # crude initial targets obtained by averaging without aligning
xs = np.arange(siz)
fall, axall = plt.subplots(1,3,figsize=(6,2),sharex=True)
for whichnoise in [0,1,2]: # the least noisy dataset is 2
    SNR = 1/noiselevels[whichnoise]**2
    print ("SNR=", SNR)
    images = samples[:,:,whichnoise]
    
    f,ax = plt.subplots(1,Nrefine,figsize=(15,2),sharex=True)
    running = avgs[:,whichnoise] # initial guess
    ax[0].plot(running)
    ax[0].set_title('SNR='+str(SNR))
    LLlist = []
    for refine in range(1,20):
        newrunning,LL = reestimate([running,noiseSDs[whichnoise],shiftSD,0], images)
        running = newrunning[0]
        print(LL)
        LLlist += [LL]
        if refine%4==1: #and refine<Nrefine*4-6: 
            npl = int((refine+3)/4)
            ax[npl].plot(running)
            ax[npl].set_title(str(refine))
    plt.figure()
    plt.plot(LLlist)
    plt.xlabel('refinement step'); plt.ylabel('LL'); 
    plt.title('SNR='+str(SNR))
    axall[whichnoise].plot(running)
    axall[whichnoise].set_title('SNR='+str(SNR))
