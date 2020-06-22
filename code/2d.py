from scipy.ndimage.interpolation import rotate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mrc
import pyfftw, numpy
import pyfftw.interfaces.numpy_fft
from numba import jit

def do_complex_rotate(arr,angle,rotate_func=rotate, **kwargs):
  r = rotate(np.real(arr),angle=angle, reshape=False, **kwargs)
  i = rotate(np.imag(arr),angle=angle, reshape=False, **kwargs)
  return(r+i*1j)

 def cmask(index,radius,array):
  a,b = index
  nx,ny = array.shape
  y,x = np.ogrid[-a:nx-a,-b:ny-b]
  mask = x*x + y*y <= radius*radius
  return(mask)

 def do_2dplot(arr):
  plt.imshow(arr, cmap='gray')

def do_1dplot(arr,idx=None,**kwargs):
  sr = pd.Series(arr)
  if idx is not None: sr.index=idx
  sr.plot(**kwargs)

def log_abs(arr):
  return(np.log(1+np.abs(arr)))

def fft2d(arr2d,mode,numpy_fft=pyfftw.interfaces.numpy_fft):
  '''
  we apply an alterating +1/-1 multiplicative before we go to/from Fourier space. 
  Later we apply this again to the transform.
  '''
  assert arr2d.ndim == 2
  n1,n2 = arr2d.shape
  assert n1==n2
  arr2d = neg_pos(arr2d.copy())
  if mode=='f':
    arr2d_f = numpy_fft.fftn(arr2d.reshape(-1,n1,n1),axes=(-2,-1))
    arr2d_f /= n1
  elif mode=='i':
    arr2d_f = numpy_fft.ifftn(arr2d.reshape(-1,n1,n1),axes=(-2,-1)).real
    arr2d_f *= n1
  
  arr2d_f = neg_pos(arr2d_f.reshape(n1,n1).copy())
  return(arr2d_f)

def do_fft(arr2d,**kwargs):
  return(fft2d(arr2d,mode='f',**kwargs))

def do_ifft(arr2d,**kwargs):
  return(fft2d(arr2d,mode='i',**kwargs))

@jit
def neg_pos(arr2d):
  '''
  each pixel switches from positive to negative in checker board pattern
  '''
  assert arr2d.ndim == 2
  for r in range(arr2d.shape[0]):
    for c in range(arr2d.shape[1]):
      if (r+c)%2:
        arr2d[r,c] *= -1
  return(arr2d)


def comp_corr(a,b):
  '''
  cross correlation in Fourier space
  the cross correlation of the real and imaginary parts and adds them up
  '''
  corr = np.multiply(np.real(a), np.real(b)).sum() + np.multiply(np.imag(a), np.imag(b)).sum()
  return(corr)