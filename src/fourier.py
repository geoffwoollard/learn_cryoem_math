import pyfftw
import numpy as np
import pyfftw.interfaces.numpy_fft
from numba import jit

@jit
def make_neg_pos_3d(arr3d):
  R1, R2, R3 = arr3d.shape
  neg_pos_3d = np.ones(arr3d.shape)
  for r1 in range(R1):
    for r2 in range(R2):
      for r3 in range(R3):
        if (r1+r2+r3)%2:
          neg_pos_3d[r1,r2,r3] *= -1.0
  return(neg_pos_3d)

def fft3d(arr3d,mode,neg_pos_3d=None,numpy_fft=pyfftw.interfaces.numpy_fft,only_real=False):
  if neg_pos_3d is None:
    neg_pos_3d = make_neg_pos_3d(arr3d)
  if arr3d.shape[0]%4 != 0:
      neg_pos_3d *= -1

  arr3d_f = numpy_fft.fftn(neg_pos_3d*arr3d)
  
  if mode == 'f':
    arr3d_f /= np.sqrt(np.prod(arr3d_f.shape))
  elif mode == 'i':
    arr3d_f *= np.sqrt(np.prod(arr3d_f.shape))

  if only_real:
    arr3d_f = arr3d_f.real

  arr3d_f *= neg_pos_3d
  return(arr3d_f)

def do_fft(arr3d,d=3,only_real=False,**kwargs):
  assert d == 3, 'only 3d implemented'
  return(fft3d(arr3d,mode='f',**kwargs))

def do_ifft(arr3d,d=3,only_real=True,**kwargs):
  assert d == 3, 'only 3d implemented'
  return(fft3d(arr3d,mode='i',**kwargs))
