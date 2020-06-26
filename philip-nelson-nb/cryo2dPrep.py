#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:32:27 2019
@author: pcn  cryo2dPrep.py create fake 2D electron-microscopy data and write file 2dImages.npz
Python 3.7
Description: output file contains images at each of 3 noise levels, all 
 jittered by the same amount in x,y and then randomly rotated about origin:
  4D array samples[which image, x, y, which_noiselevel]
  array noiselevels[which_noiselevel] = SD of noise used to corrupt samples
  value shiftSD = jitter in pixels used to corrupt samples

creates images
"""
import numpy as np; import matplotlib.pyplot as plt
from numpy.random import randn, random
from scipy.ndimage.interpolation import rotate
from scipy.interpolate import griddata

plt.close('all')
siz = 85    # size of each 2D "image"
if siz%2==0: raise ValueError('siz should be odd')
Nsamp = 1500 # how many samples taken; =1500 in notes
shiftSD = 6 # how many pixels to jitter the fake data, normally 6
do_rotations = True  # False if you only want shifts (create image 2dNoisyJit)
    # True if you want shifts and rotations (create image 2dNoisyJitRot)
myImage = 'ChineseLifeTilt85.tif' # alts 'ChineseLifeTilt85.tif' 'ChineseLife85.tif'; 'test2D.tif' 

noiselevels = np.array([6, 3, 1]) # = 1/sqrt(SNR) because signal normalized to var=1; =[6, 3, 1] in notes
print(myImage, noiselevels, Nsamp, do_rotations)

def myimshow(y,ax):
    ax.imshow(y,'gray')
    ax.axis('off')

# get the "true image"
ys = plt.imread(myImage).astype('float')
print(ys.shape)
f,ax = plt.subplots(1,1,figsize=(1.5,1.5))
signalsig = np.var(ys)
ys = ys/np.sqrt(signalsig) # normalize to unit var
print('check=1: ', np.var(ys))
myimshow(ys,ax)
ys2 = ys.reshape((1,siz,siz))

noises = randn(Nsamp, siz, siz)
samples = np.zeros((Nsamp, siz, siz, 3))


# display a few examples of noisy data and show noise suppression by avg:
fig,ax = plt.subplots(3,4,figsize=(8,7),sharex=True, sharey=True)
for k in range(len(noiselevels)):
    noi = noiselevels[k]
    samples[:,:,:,k] = (ys2 + noi*noises) # sd of the noise is noiselevels[k]. broadcast ys.
    ax[k,0].set_title('SNR='+str(1/noiselevels[k]**2)+' unjittered, 3 exemplars')
    for which in [0,1,2]:
        myimshow(samples[which,:,:,k], ax[k,which])
    myimshow(np.mean(samples[:,:,:,k],axis=0), ax[k,3])
    ax[k,3].set_title('mean')


# now jitter the noisy data in x,y
jitters = randn(Nsamp,2)*shiftSD # translations
rotters = random(Nsamp)*2*np.pi  # rotations
buffers = int(2*(np.sqrt(2)-1)) + 4*shiftSD   # padding to make sure rotations fit
padsiz = siz + 2*buffers
xs = np.arange(siz)
xvals, yvals = np.meshgrid(xs, xs, indexing='ij')
padxs = np.arange(padsiz) - buffers
padxvals, padyvals = np.meshgrid(padxs, padxs, indexing='ij')
lpadxs = padxvals.ravel(); lpadys = padyvals.ravel() # linear listings
half = (siz-1)/2
paddedx = np.concatenate((np.zeros((buffers,siz)), ys, \
                         np.zeros((buffers,siz))), axis=0)
padded = np.concatenate((np.zeros((padsiz,buffers)), paddedx, \
                         np.zeros((padsiz,buffers))), axis=1)
for k in range(3):  # which noise level
    print('noise',k)
    for j in range(Nsamp):  # which image
        if j%100==0: print(j)
        angle = rotters[j]  # how much to rotate this instance
        if do_rotations:
            rotated = rotate(padded, angle*360/(2*np.pi), reshape=False)
        else:
            rotated = padded
        lpadded = rotated.ravel() # linear listing
        jitxvals = xvals + jitters[j,0]
        jityvals = yvals + jitters[j,1]
        samples[j,:,:,k] = noises[j,:,:]*noiselevels[k] + \
         griddata((lpadxs, lpadys), lpadded, (jitxvals.ravel(), jityvals.ravel()), \
         method='cubic', fill_value=0., rescale=False).reshape((siz,siz))
            

fig,ax = plt.subplots(3,4,figsize=(8,7), sharex=True, sharey=True)
for k in range(len(noiselevels)):
    noi = noiselevels[k]
    ax[k,0].set_title('SNR='+str(1/noiselevels[k]**2)+' jittered, 3 exemplars')
    for which in [0,1,2]:
        myimshow(samples[which,:,:,k], ax[k,which])
    myimshow(np.mean(samples[:,:,:,k],axis=0), ax[k,3])
    ax[k,3].set_title('mean')

np.savez('2dImages',samples=samples,noiselevels=noiselevels,shiftSD=shiftSD,myImage=myImage)