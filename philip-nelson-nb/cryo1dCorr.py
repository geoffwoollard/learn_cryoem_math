#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 07:52:33 2019
@author: pcn  cryo1dCorr.py attempt to align noisy samples by cross-correlation to their average
Python 3.7
Description: 
"""
import numpy as np; import matplotlib.pyplot as plt
from scipy import signal

# parameters
iexamine = 1 # choose one image to show details about

plt.close('all')
a=np.load('1dImages.npz')
shiftSD = a['shiftSD'] # amplitude of jitter in the fake data
samples=a['samples']
noiselevels=a['noiselevels']
Nsamp,siz,Nnois = samples.shape
halfsize = int((siz-1)/2)

avgs = np.mean(samples, axis=0) # initial targets
xs = np.arange(siz)

# test to make sure we understand the options for corr
test=avgs[:,2]
test2 = np.append(test[2:],[0.,0.]) #second has been shifted left 2 slots
mycor = signal.correlate(test,test2,mode='same')
themax = np.argmax(mycor)
plt.plot(mycor); plt.title("Same:"+str(len(mycor))+": "+str(themax))
# preceding exercise shows that corr' peaks at (M-1)/2+(left shift of 2nd series wrt 1st)
mycor = signal.correlate(test,test2,mode='full')
themax = np.argmax(mycor)
plt.figure(); plt.plot(mycor); plt.title("Full:"+str(len(mycor))+": "+str(themax))
# preceding exercise shows that full corr peaks at M-1+(left shift of 2nd series wrt 1st)

# align by shifting test left by 2 slots:
leftshift = themax - (siz-1)
corrected = np.append(test[leftshift:],np.zeros(leftshift))
plt.figure(); plt.plot(xs,test2,xs,corrected+.02)

fig,ax = plt.subplots(3,4,figsize=(8,7))
for k,noi in enumerate(noiselevels):
    newavg = np.zeros(siz) # accumulate average of shifted samples
    thisavg = avgs[:,k]   #  try to align to this target
    for i in range(Nsamp):
        s = samples[i,:,k] 
        c = signal.correlate(s, thisavg, mode='same')
        if i==0:
            ax[k,0].set_ylabel('SNR='+str(1/noiselevels[k]**2))
            ax[k,0].set_title('target')
            ax[k,0].plot(thisavg)
            ax[k,1].set_title('simulated data')
            ax[k,1].plot(s)
            ax[k,2].set_title('correlation')
            ax[k,2].plot(c)
        leftshift = np.argmax(c) - halfsize # s should be left-shifted by this much to align with thisavg
        if i==iexamine: 
            print('image number ',i,'; best leftshift=', leftshift)
        if leftshift>0: #shift left
            s = np.append(s[leftshift:],np.zeros(leftshift))
        elif leftshift<0: #shift right
            rsh = -leftshift
            s=np.append(np.zeros(rsh),s[:-rsh])
#        ax[3].plot(xs,thisavg,xs,s); ax[3].set_title('shifted') # temp
        newavg += s/Nsamp # accumulate average of shifted samples
    ax[k,3].plot(newavg)
    ax[k,3].set_title('new average')
