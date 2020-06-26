#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on 19 Aug 07:52:33 2019
@author: pcn  cryo2dCorr.py attempt to align noisy samples by cross-correlation to 
a starting image. 
Description: Although each sample is unblurred, each is corrupted by additive gaussian noise and 
unknown random shifts and rotations.
Python 3.7
"""
import numpy as np; import matplotlib.pyplot as plt
from scipy import signal
from scipy.ndimage.interpolation import rotate


# Parameters:
angleres = 180 # angular resolution for our correlation search (divide circle into this many divisions)
limit_samp = 10000 # restrict to this many samples for debugging, or all samples in the file, whichever is less
mytarget = 'test2Dblur.tif' #'ChineseLifeTiltBlur3px.tif' # which image to use as the template
# alts 'test2Dblur.tif'; 'ChineseLifeBlur2px.tif' 'ChineseLifeBlur3px.tif' 'ChineseLifeTiltBlur3px.tif'

# get the "starting image" (wish to improve it using the "experimental" data):
target = plt.imread(mytarget).astype('float')

# Get fake "experimental" data:
plt.close('all')
a=np.load('2dPrepChinese6-3-1-1500NewRotation.npz')# ('2dImages.npz')
shiftSD = a['shiftSD']; print('jitter=',shiftSD)  # amplitude of jitter in the fake data
myImage = a['myImage']; print(myImage, mytarget)
samples=a['samples'] # this is the simulated data
noiselevels=a['noiselevels']; print('noiselevels=',noiselevels,'; SNRs=', 1/noiselevels**2)
Nsamp,siz,_,Nnois = samples.shape
if Nsamp>limit_samp:  # for debugging
    samples= samples[:limit_samp,:,:,:]
    Nsamp = limit_samp
halfsize = int((siz-1)/2)

print('target image: ', mytarget, '; ', Nsamp, ' samples')

#%%
def myimshow(y,ax):
    ax.imshow(y,'gray')
    ax.axis('off')
#%%
# Make a table of all rotated views of the target
targetsrot = np.zeros((siz,siz,angleres)) # preallocate
rotters = np.linspace(0,2*np.pi,angleres+1)[:-1] # rotations
xs = np.arange(siz)
half = siz/2  # not integer but that's ok

xvals, yvals = np.meshgrid(xs,xs, indexing='xy')
lpadxs = xvals.ravel(); lpadys = yvals.ravel() # linear listings

#tmp = interp2d(padxs-buffers, padxs-buffers, padded, kind='cubic')
fig,ax = plt.subplots(1,5,figsize=(7,2))
for m,angle in enumerate(rotters):
    targetsrot[:,:,m] = \
             rotate(target, -angle*360/(2*np.pi), reshape=False)
    if m<5: myimshow(targetsrot[:,:,m], ax[m])

#%%
iexamine = 1 # which image to display

fig,ax = plt.subplots(3,6,figsize=(10,7))
for k,noi in enumerate(noiselevels):
    print('noise level=',k)
    newavg = np.zeros((siz,siz)) # accumulate average of aligned samples
    for i in range(Nsamp):
#        leftshifts = np.zeros(angleres) # allocate for best shifts for each view
#        upshifts = np.zeros(angleres); 
        s = samples[i,:,:,k]
        bestcorvalue = -1e20 # minus infinity
        for m in range(angleres): # check all orientations
            thistarget = targetsrot[:,:,m]
            c = signal.correlate(s, thistarget, mode='same')
            if (i==iexamine) and (m==0):
                ax[k,0].set_title('target')
#                ax[k,0].set_ylabel('SNR = '+str(1/thisnoise**2))
                myimshow(targetsrot[:,:,0], ax[k,0])
                ax[k,1].set_title('sim. data,SNR='+str(1/noi**2))
                myimshow(s, ax[k,1])
                ax[k,2].set_title('corr,rand.orient.')
                myimshow(c, ax[k,2])
                print('best cor, random orient.=',np.amax(c))
            shifmax = np.unravel_index(np.argmax(c, axis=None), c.shape) # get row and col of best shift
            thisbestcor = np.amax(c)
            if thisbestcor > bestcorvalue:   # best found so far
                bestcors = c # later display best one
                bestcorvalue = thisbestcor
                bestrot = rotters[m]
                bestupshift = shifmax[0] - halfsize # thistarget is upward-shifted by this much w.r.t. s
                # hence s should be upward-shifted by this much to align with thistarget
                bestleftshift = shifmax[1] - halfsize # thistarget is left-shifted by this much w.r.t. s
                # hence s should be left-shifted by this much to align with thistarget
        if i==iexamine: 
            print('image number ',i,' best orient m=',bestrot,'; best upshift=',bestupshift, \
                  '; best leftshift=', bestleftshift)
            print('max cor, best orient=', bestcorvalue)
            myimshow(bestcors , ax[k,3])
            ax[k,3].set_title('corr,best orient.')
        #we now know the best-matching orientation and shift, so align image i:
        if bestleftshift>0: #shift image left
            s = np.hstack(( s[:,bestleftshift:], np.zeros((siz,bestleftshift)) ))
        elif bestleftshift<0: #shift right
            rsh = -bestleftshift
            s = np.hstack((np.zeros((siz,rsh)), s[:,:-rsh]))
            
        if bestupshift>0: #shift up
            s = np.vstack(( s[bestupshift:,:], np.zeros((bestupshift,siz)) ))
        elif bestupshift<0: #shift down
            dsh = -bestupshift
            s = np.vstack((np.zeros((dsh,siz)), s[:-dsh,:]))
        s = rotate(s, bestrot*360/(2*np.pi), reshape=False)
        if i==iexamine:
            myimshow(s , ax[k,4])
            ax[k,4].set_title('best align')
            
        newavg += s/Nsamp # accumulate average of shifted samples
    myimshow(newavg, ax[k,5])
    ax[k,5].set_title('new average')
