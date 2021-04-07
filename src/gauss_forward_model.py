import numpy as np
from numba import njit, prange

def make_gauss_2d(xv,yv,mu,sigma):
  g = np.exp(-( (xv-mu[0])**2 + (yv-mu[1])**2)  /(2*sigma**2) )
  return(g)

def make_map_3d(atoms,xyz,sigma):
  C = 1/np.sqrt(2*np.pi*sigma**2)
  N = atoms.shape[1]
  diff = xyz.reshape(-1,3,1) - atoms[:3,:].reshape(1,3,-1)
  a = -1/(2*sigma**2)
  map_3d = (C**3*np.exp(a*((diff**2).sum(1))).sum(1).reshape(N,N,N))
  return(map_3d)

@njit(parallel=True)
def parallel_add_patch_including_diff(xy,N,atoms,idx,n_trunc,sigma):
  nt_ = (n_trunc-1)//2
  a = -1/(2*sigma**2)
  g_2d = np.zeros(N*N)
  for i in prange(idx.shape[0]): # loop over atoms
    patch_line = np.arange(idx[i]-nt_,idx[i]+nt_+1,1)
    for y_line in prange(-nt_,+nt_+1,1):
      one_gauss_patch_y_line_idxs = patch_line + y_line*N

      diffx = xy[one_gauss_patch_y_line_idxs,0] - atoms[0,i]
      diffy = xy[one_gauss_patch_y_line_idxs[0],1] - atoms[1,i] # all ys the same
      d2i = diffx*diffx+diffy*diffy
      gi_y_line = np.exp(a*d2i)
      g_2d[one_gauss_patch_y_line_idxs] += gi_y_line

  return(g_2d)

@njit(parallel=True)
def parallel_add_patch_from_idx(idx,d2,N,n_trunc,sigma):
  # https://numba.pydata.org/numba-doc/latest/user/parallel.html njit and prange with +=
  nt_ = (n_trunc-1)//2
  a = -1/(2*sigma**2)
  g_2d = np.zeros(N*N)
  for i in prange(idx.shape[0]): # loop over atoms
    patch_line = np.arange(idx[i]-nt_,idx[i]+nt_+1,1)
    for y_line in prange(-nt_,+nt_+1,1):
      one_gauss_patch_y_line_idxs = patch_line + y_line*N
      gi_y_line = np.exp(a*d2[one_gauss_patch_y_line_idxs,i])
      g_2d[one_gauss_patch_y_line_idxs] += gi_y_line

  return(g_2d)

def make_proj_mask(atoms, xy, sigma, n_trunc,parallel_diff=True):
  N = np.sqrt(xy.shape[0]).astype(int)
  X = np.round(atoms[0]).astype(np.int32) + N//2
  Y = np.round(atoms[1]).astype(np.int32) + N//2
  idx = X+N*Y
  if parallel_diff:
    g_2d = parallel_add_patch_including_diff(xy,N,atoms,idx,n_trunc,sigma).reshape(N,N)
  else:
    diff = xy.reshape(-1,2,1) - atoms[:2,:].reshape(1,2,-1) # proj indep of z
    d2 = (diff**2).sum(axis=1)
    g_2d = parallel_add_patch_from_idx(idx,d2,N,n_trunc,sigma).reshape(N,N)
  return(g_2d)