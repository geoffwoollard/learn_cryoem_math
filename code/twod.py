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


def comp_corr(a,b,axis):
  '''
  cross correlation in Fourier space
  the cross correlation of the real and imaginary parts and adds them up
  '''
  corr = np.multiply(np.real(a), np.real(b)).sum(axis=axis) + np.multiply(np.imag(a), np.imag(b)).sum(axis=axis)
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

def do_2d_align_poisson(X,
  n_A_updates,
  A_prev=None,
  noise_param_d=None, #lam_k, sigma
  deg_step=None,
  shift_span=0,
  sigma_shift=np.inf,
  bool_circle_mask=None,
  do_plot=True,
  figsize=(16,32),
  do_log=False,
  X_aligned=None,
  A=None,
  stats='poisson'):

  if A_prev is None: A_prev = X.mean(0)
  assert A_prev.shape == X.mean(0).shape
  A_next = A_prev.copy()

  small_N = X.shape[0]
  nx = X.shape[-1]

  # mask
  if bool_circle_mask is None:
    bool_circle_mask = ~cmask(index=(nx//2,nx//2),radius=int(nx//2*0.85),array=np.ones_like(A_prev)).astype(np.bool)

  # shifts
  shifts_r = np.arange(-shift_span,shift_span+1, dtype=np.int32)
  shifts_c = np.arange(-shift_span,shift_span+1, dtype=np.int32)
  if shift_span==0: 
    shifts_r = np.array([0])
    shifts_c = np.array([0])

  # angles
  angles = np.arange(1,360,deg_step)
  if deg_step is None:
    angles = np.array([0])

  # initialize
  A_align = np.zeros((nx,nx,angles.shape[0],shifts_r.shape[0],shifts_r.shape[0] ))
  x_aligned = np.zeros_like(A_align)
  best_X = np.zeros((n_A_updates,small_N,nx,nx))
  best_angles, best_shift_rs, best_shift_cs = np.zeros((n_A_updates,small_N)), np.zeros((n_A_updates,small_N)), np.zeros((n_A_updates,small_N))
  LL = np.zeros((n_A_updates,small_N))
  A_nexts = np.zeros((n_A_updates,) + A_next.shape)

  if stats == 'gaussian':
      sigma = noise_param_d['sigma']

  extra_plot_n=4
  for obj in [A,X_aligned]:
    if obj is None: plot_n -= 1

  if do_plot: fig, axes = plt.subplots(min(10,small_N)+extra_plot_n, n_A_updates,figsize=figsize)

  for c in range(n_A_updates):
    if do_log: print(c)
    
    ll=0
    A_prev = A_next.copy()
    A_next = np.zeros_like(A_prev)

    if do_plot: axes[0,c].imshow(A_prev,cmap='gray') ; axes[0,c].set_axis_off()

    # reference alignments of template
    for shift_r_idx, shift_r in enumerate(shifts_r):
      A_shift_r = shift_zeropad_axis(A_prev,shift=shift_r,axis=0)
      for shift_c_idx, shift_c in enumerate(shifts_c):
        A_shift_r_c = shift_zeropad_axis(A_shift_r,shift=shift_c,axis=1)
        for angle_idx, angle in enumerate(angles):
          A_shift_r_c[bool_circle_mask] = 0 # TODO test if interpolation different with windowing
          A_align[:,:,angle_idx,shift_r_idx,shift_c_idx] = rotate(A_shift_r_c,angle=angle, reshape=False)
    
    # terms that only depends on A
    if stats=='poisson':
      lam_k = noise_param_d['lam_k']
      negs = A_align[~bool_circle_mask][A_align[~bool_circle_mask] < 0]
      if negs.size < 0:
        negs = A_align[~bool_circle_mask][A_align[~bool_circle_mask] > 0].min() # hack to clip to smallest non zero value
      log_lam = np.log(A_align[~bool_circle_mask]+lam_k)
      # table of norms
      log_etolam = -(A_align[~bool_circle_mask].sum(axis=0)+(lam_k)*A_align[~bool_circle_mask].size) # the mask collapses the two xy image axes into one
    elif stats == 'gaussian':
      A_aligned_norm = np.linalg.norm(A_align[~bool_circle_mask],axis=0)
      A_aligned_norm_  = -(2*sigma**2)**-1*A_aligned_norm
    else:
      assert False, 'only poisson and gaussian stats implemented'

    # pdf shifts, shift prior
    log_prior_shift = np.zeros_like(A_align[0,0])
    for shift_r_idx in range(shifts_r.shape[0]):
      for shift_c_idx in range(shifts_c.shape[0]):
        q2 = shifts_r[shift_r_idx]**2+shifts_c[shift_c_idx]**2
        log_prior_shift[:,shift_r_idx,shift_c_idx] = -q2/(2*sigma_shift**2)

    r=0
    sigma_update = 0
    for i in range(small_N):
      #print('image %i'%i)
      x = X[i]
          
      #Ki, gi
      if stats == 'poisson':
        log_lamtok = x[~bool_circle_mask].reshape(x[~bool_circle_mask].shape+(1,1,1,))*log_lam
        log_gi_align = log_lamtok.sum(axis=0) + log_etolam + log_prior_shift

      elif stats == 'gaussian':
        for angle_idx, angle in enumerate(angles):
          newshape = x.shape + tuple(np.ones(A_align.ndim-2,dtype=np.int32))
          corr_A_x[angle_idx] = comp_corr(A_align[:,:,angle_idx][~bool_circle_mask],
                                        x.reshape(newshape)[~bool_circle_mask],
                                        axis=0) # vectorized over alignment, axis 0 is pixels (flattened from bool_circle_mask)
        corr_A_x_ = sigma**-2*corr_A_x
        log_gi_align = A_aligned_norm_ + corr_A_x_ + log_prior_shift
      else:
        assert False, 'only poisson and gaussian stats implemented'
      
      Ki = log_gi_align.max()
      log_gi_align_stable = log_gi_align - Ki
      gi_stable = np.exp(log_gi_align_stable, dtype=np.float128)
     
      # Ui
      gisum = gi_stable.sum()
      if not np.isclose(gisum, 0): 
        Ui = gisum**-1
        # log lik
        if np.isfinite(sigma_shift):
          log_sigma_shift = np.log(sigma_shift)
        else:
          log_sigma_shift = 0

        ll += -np.log(Ui) + Ki - sum_ln_factorial(x) -0.5*np.log(2*np.pi) - log_sigma_shift
        
      else: 
        Ui=0

      LL[c,i] = np.log(-ll)

      # update noise model
      if stats == 'gaussian':
        #newshape = x.shape + tuple(np.ones(A_align.ndim-2,dtype=np.int32))
        diff  = np.subtract(A_align,x.reshape(newshape))
        errors = np.linalg.norm(diff[~bool_circle_mask],axis=0)**2
        sigma2_i = (gi*errors).sum()
        sigma2_i /= gi.sum()
        sigma_update += np.sqrt(sigma2_i)
        if do_log: print('sigma_i',np.sqrt(sigma2_i))

      # rev alignment
      x_aligned = comp_x_aligned(x,A_align,angles,shifts_r,shifts_c)

      # point estimate of best angle
      angle_idx_best, shift_r_idx_best,shift_c_idx_best = np.unravel_index(np.argmax(gi_stable, axis=None), gi_stable.shape)
      best_angles[c,i] = angles[angle_idx_best]
      best_shift_rs[c,i] = shifts_r[shift_r_idx_best]
      best_shift_cs[c,i] = shifts_c[shift_c_idx_best]
      best_X[c,i,:,:] = x_aligned[:,:,angle_idx_best,shift_r_idx_best,shift_c_idx_best]

      # Maximization (update A)
      A_next += Ui*np.multiply(gi_stable.reshape((1,1,)+gi_stable.shape),x_aligned).sum(axis=(-1,-2,-3))

      if i % np.ceil(X[:small_N].shape[0]/10) == 0: 
        if do_log: print('i = %i, ll = %.2f, A_next min=%.2f, max=%.2f' % (i,ll,A_next[~bool_circle_mask].min(),A_next[~bool_circle_mask].max()))
        if do_plot: 
          axes[r+1,c].imshow(A_next,cmap='gray')
          axes[r+1,c].set_axis_off()
          r+=1

    A_next /= small_N
    A_nexts[c] = A_next
    if stats == 'gaussian': 
      if do_log: print('sigma_update',np.sqrt(sigma_update))
      sigma = sigma_update

    if do_plot: 
      axes[r+1,c].imshow(X[:small_N].mean(0),cmap='gray') ; axes[r+1,c].set_axis_off()
      if X_aligned is not None: axes[r+2,c].imshow(X_aligned[:small_N].mean(0),cmap='gray') ; axes[r+2,c].set_axis_off()
      if A is not None: axes[r+3,c].imshow(A,cmap='gray') ; axes[r+3,c].set_axis_off()

  return(A_nexts)