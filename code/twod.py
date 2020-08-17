from scipy.ndimage.interpolation import rotate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mrc
import pyfftw, numpy
import pyfftw.interfaces.numpy_fft
from numba import jit

def shift_zeropad_axis(x,shift,axis):
  assert axis in [0,1]
  if axis == 0:
    x = np.roll(x,shift,axis=axis)
    
    if shift > 0:
      x[:shift,:] = 0
    elif shift < 0:
      x[shift:,:] = 0
    else: pass
  
  else:
    x = np.roll(x,shift,axis=axis)
    if shift > 0:
      x[:,:shift] = 0
    elif shift < 0:
      x[:,shift:] = 0
    else: pass
  return(x)

def comp_x_aligned(x,A_rot_shifted,angles,shifts_r,shifts_c):
  '''
  TODO: rewrite without A_rot_shifted, since just using for shape
  '''
  x_aligned = np.zeros(A_rot_shifted[:,:,:,:,:].shape)
  for angle_idx in range(angles.shape[0]):
    x_rot = rotate(x,angle=-angles[angle_idx],reshape=False) 
    for shift_r_idx in range(shifts_r.shape[0]):
      x_rot_shift = shift_zeropad_axis(x_rot,shift=-shifts_r[shift_r_idx],axis=0)
      for shift_c_idx in range(shifts_c.shape[0]):
          x_aligned[:,:,angle_idx,shift_r_idx,shift_c_idx] = shift_zeropad_axis(x_rot_shift,shift=-shifts_c[shift_c_idx],axis=1)
  return(x_aligned)

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

def fft2d(arr2d,mode,numpy_fft=pyfftw.interfaces.numpy_fft,only_real=True):
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
    if only_real:
      arr2d_f = numpy_fft.ifftn(arr2d.reshape(-1,n1,n1),axes=(-2,-1)).real
    else:
      arr2d_f = numpy_fft.ifftn(arr2d.reshape(-1,n1,n1),axes=(-2,-1))
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

def simulate_data(image_2d,psize_A,N_particles=500,df_low=1e4,df_high=2e4,snr=3,bool_circle_mask=None):
  assert image_2d.ndim == 2
  nx = image_2d.shape[0]
  image_2d_f = do_fft(image_2d) # ground truth image

  true_angles = np.random.uniform(low=0,high=360, size=N_particles)

  ctf_2ds = np.zeros((N_particles,nx,nx))
  dfs = np.random.uniform(low=df_low,high=df_high,size=N_particles)
  s, a = ctf.ctf_freqs(image_2d.shape,d=1/psize_A)

  noise_std = image_2d_f.std() / snr

  images_observed = np.zeros((N_particles,nx,nx), dtype=np.complex64)

  sim_params_d = defaultdict(list)

  for i in range(N_particles):
    # rotate
    image_2d_f_rot = do_complex_rotate(image_2d_f,angle=true_angles[i])
    if bool_circle_mask is not None:
      image_2d_f_rot[bool_circle_mask]=0
    # ctf
    ctf_2ds[i,:,:] = np.fft.fftshift(ctf.eval_ctf(s, a, def1=dfs[i], def2=dfs[i], angast=0, phase=0, kv=300, ac=0.1, cs=2.0, bf=0, lp=0))
    image_2d_f_rot_ctf = image_2d_f_rot*ctf_2ds[i]
    # noise
    noise = np.random.normal(loc=0,scale=noise_std,size=nx*nx).reshape(nx,nx)
    images_observed[i,:,:] = image_2d_f_rot_ctf + noise

  sim_params_df = pd.DataFrame({'df1':dfs,'df2':dfs, 'pose2D':true_angles})
  return(images_observed,sim_params_df)

def sum_ln_factorial(x):
  '''
  sum log of lectron count factorial over pixels
  precompute pixel_values and their counts 
  sum_a ln Xia!
  '''
  
  value_counts = pd.Series(x.astype(np.uint64)[x>1]).value_counts() # 0 and 1 contribute nothing to sum
  pixel_values = value_counts.index.values.astype(np.uint64)
  counts = value_counts.values.astype(np.uint64)
  lnxia = 0
  for pixel_value, count in zip(pixel_values, counts):
    lnxia += count*np.math.factorial(pixel_value)
  return(lnxia)

@jit
def rotate_bi(arr,angle):
  '''
  rotation with bilinear interpolation
  cartesian to cartesian, no need to interpolate whole image onto polar. when convert to polar, just need to add angle of rotation
  see http://polymathprogrammer.com/2008/10/06/image-rotation-with-bilinear-interpolation/
  '''
  arr_ = np.zeros_like(arr)
  I,J = arr_.shape # loop over destination
  
  # rotation angle (rad)
  rad = angle*np.pi/180

  for i in range(I):
    y = I//2 - i
    for j in range(J):
      # raster to cartesian
      x = j - J//2
      

      # cartesian to polar
      r = np.sqrt(x*x+y*y)
      t = np.arctan2(y,x)


      #polar to cartesian (in new ref frame)
      x_ = r*np.cos(t-rad)
      y_ = r*np.sin(t-rad)

      # cartesian to raster
      j_ = x_ + J//2
      i_ = I//2 - y_

      # floor and ceil (of one pixel)
      i_floor = int(np.floor(i_))
      j_floor = int(np.floor(j_))
      i_ceil = int(np.ceil(i_))
      j_ceil = int(np.ceil(j_))

      # check bounds of the pixel are within image
      test = (i_floor < 0 or j_floor < 0 or i_ceil >= I or j_ceil >= J or i_floor >= I or j_floor >= J or i_ceil < 0 or j_ceil < 0)
      if test: continue
    
      di = i_ - i_floor
      dj = j_ - j_floor

      # corners (of one pixel)
      top_left = arr[i_floor,j_floor]
      top_right = arr[i_floor,j_ceil]
      bot_left = arr[i_ceil,j_floor]
      bot_right = arr[i_ceil,j_ceil]


      # linearly interpolate horizontally between top neighbours
      top = (1-dj)*top_left + dj*top_right

      # linearly interpolate horizontally between bottom neighbours
      bot = (1-dj)*bot_left + dj*bot_right

      # linearly interpolate vertically between top and bottom interpolated results
      p_ = (1-di)*top + di*bot

      arr_[i,j] = p_ # this doesn't need to be rounded since not 8 bit
  return(arr_)