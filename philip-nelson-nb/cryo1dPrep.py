#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:32:27 2019
@author: pcn  cryo1dPrep.py create fake 1D electron-microscopy data and write file 1dImages.npz
Python 3.7
Description: output file contains images at each of 3 noise levels, all 
 jittered by the same amount in x,y:
array samples[which_image, x, which_noiselevel]
array noiselevels[which_noiselevel] used to corrupt samples
value shiftSD = jitter in pixels used to corrupt samples
"""
import numpy as np; import matplotlib.pyplot as plt
from numpy.random import randn
from scipy.interpolate import interp1d

plt.close('all')

# parameters
iexamine = 1 # choose one image to show details about
siz = 85    # size of each 1D "image"
if siz%2==0: raise ValueError('siz should be odd')
Nsamp = 1500 # how many samples taken; later much bigger
shiftSD = 6 # how many pixels to jitter the fake data

noiselevels = np.array([2,1,.25]) # = 1/sqrt(SNR) because signal normalized to var=1; [2,1,.25] in notes 

# make the "true image": square peak, gaussian peak, triangular peak
xs = np.arange(siz)             # pixel addresses
ys = np.exp(-(xs-siz/2)**2/23)  # true image
ys[10:15] = .8
peak = np.linspace(0,5,5)/6
ys[60:65]=peak
ys[64:68]=peak[-1:0:-1]
signalsig = np.var(ys)
ys = ys/np.sqrt(signalsig) # normalize to unit var
print('check=1: ', np.var(ys))
plt.plot(xs,ys,'.--')
plt.xlabel(r'$x,\ \mathrm{pixels}$')
plt.ylabel(r'$\mathrm{intensity\ \ [}\mathsf{a.u.}\mathrm ]$')
plt.figure() # now represent as a raster:
plt.imshow(ys*np.ones((10,1)), cmap='gray')
noises = randn(Nsamp, siz)
samples = np.zeros((Nsamp, siz, 3))

# display a few examples of noisy data and show noise suppression by avg:
fig,ax = plt.subplots(3,4,figsize=(8,7),sharex=True)
for k in range(len(noiselevels)):
    noi = noiselevels[k]
    samples[:,:,k] = (ys + noi*noises) # sd of the noise is noiselevels[k]
    ax[k,0].set_ylabel('SNR='+str(1/noiselevels[k]**2))
#    plt.title('noise='+str(noi)+' unjittered')
    for which in [0,1,2]:
        ax[k,which].plot(xs,samples[which,:,k])
    ax[k,3].plot(xs, np.mean(samples[:,:,k],axis=0))

# now jitter the noisy data by the same random displacements for each noise level
jitters = randn(Nsamp)*shiftSD
print('shifting sample ',iexamine,' left by ', jitters[iexamine])
padded = np.hstack((np.zeros(4*shiftSD), ys, np.zeros(4*shiftSD)))
tmp = interp1d(np.arange(siz+8*shiftSD)-4*shiftSD, padded, kind='cubic')

for k in range(3):
    for j in range(Nsamp):
        samples[j,:,k] = tmp(xs+jitters[j]) + noises[j,:]*noiselevels[k] # shift sample j left by jitters[j]

fig,ax = plt.subplots(3,4,figsize=(8,7),sharex=True)
for k in range(len(noiselevels)):
    noi = noiselevels[k]
    ax[k,0].set_ylabel('SNR='+str(1/noiselevels[k]**2))
    plt.title('noise='+str(noi)+' jittered')
    for which in [0,1,2]:
        ax[k,which].plot(xs,samples[which,:,k])
    ax[k,3].plot(xs, np.mean(samples[:,:,k],axis=0))

np.savez('1dImages',samples=samples,noiselevels=noiselevels,shiftSD=shiftSD)