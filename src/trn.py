import numpy as np
import coords

def trn_rm0(map_3d,M,random_seed=None):
  """
  Parameters
  ----------
    map_3d : numpy.ndarray, shape (N,N,N)
      3d map
    M : int
      number of pseudo atoms
      random_seed : int

  Returns
  ------- 
    rm0 : numpy.ndarray, shape (M,3)
      initial sample of M points from map_3d. probability proportional to density
    map_3d_flat : numpy.ndarray, shape (N**3,)
      flattened map_3d, normalized to sum to 1
    map_3d_idx : numpy.ndarray, shape (N**3,)
      array of integers from 0 to N**3-1
    xyz : numpy.ndarray, shape (M,3)
      coordinates of voxels in map_3d

  """
  assert np.unique(map_3d.shape).size == 1, 'map must be cube, not non-cubic rectangular parallelepiped'
  N = map_3d.shape[0]
  assert N%2 == 0, 'N must be even'
  map_3d /= map_3d.sum() # 3d map to probability density
  map_3d_flat = map_3d.flatten()
  map_3d_idx = np.arange(map_3d_flat.shape[0])
  if random_seed is not None: np.random.seed(seed=random_seed)

  # this scales with M (the number of chosen items), not map_3d_idx (the possibilities to choose from)
  samples_idx = np.random.choice(map_3d_idx,size=M,replace=True,p=map_3d_flat) # chosen voxel indeces
  coords_1d = np.arange(-N//2,N//2)
  xyz = coords.coords_n_by_d(coords_1d,d=3)
  rm0 = xyz[samples_idx] # pick out coordinates that where chosen. note that this assumes map_3d_idx matches with rows of xyz
  return rm0,map_3d_flat,map_3d_idx,xyz,coords_1d

def trn_iterate(rm0,
  map_3d_flat,
  map_3d_idx,
  xyz,
  n_save=10,
  e0=0.3,
  ef=0.05,
  l0=None,
  lf=0.5,
  tf=None,
  do_log=True,
  log_n=10,
):
  """
  topology representing network
  diffuses pseudo atoms to cover probability distribution uniformly
  evolves by equation rm[t+1] = rm[t] + e exp[-km/l] * (r - rm[t]), where
    m indexes the pseudo atoms
    rm is the 3D cartesian location of the mth pseudo atom
    t is time
    r is a randomly sampled point from the initial 3d density (independent of t and therefore can be precomputed)
    e is a time depedent step size and evolves as e0(ef/e0)**(t/tf)
    km is the rank of the rm to r. how close rm is to r, compared with the other rm's, expressed as a rank
    l is a time depedent scaling of km size and evolves as l0(lf/l0)**(t/tf)
    tf is the total time steps
  see ref 
    Zhang, Y., Krieger, J., Mikulska-Ruminska, K., Kaynak, B., Sorzano, C. O. S., Carazo, J. M., … Bahar, I. (2021). 
    State-dependent sequential allostery exhibited by chaperonin TRiC/CCT revealed by network analysis of Cryo-EM maps. 
    Progress in Biophysics and Molecular Biology, 160, 104–120. 
    http://doi.org/10.1016/j.pbiomolbio.2020.08.006
  and python code in prody 
    http://prody.csb.pitt.edu/manual/reference/proteins/emdfile.html#TRNET
    https://github.com/prody/ProDy/blob/697220825ebc7498d64f4e82f53bb7ff6d98027c/prody/proteins/emdfile.py#L466

  Parameters
  ----------
    rm0 : numpy.ndarray, shape (M,3)
      initial sample of M points from map_3d. probability proportional to density
    map_3d_flat : numpy.ndarray, shape (N**3,)
      flattened map_3d, normalized to sum to 1
    map_3d_idx : numpy.ndarray, shape (N**3,)
      array of integers from 0 to N**3-1
    xyz : numpy.ndarray, shape (M,3)
      coordinates of voxels in map_3d
    n_save : int
      number of steps to save (evenly spaced out), in addition to initial step 0
    e0 : float
      initial step size
    ef : float
      final step size
    l0 : float
      initial scaling of rank. larger tightens things up (pulls together). smaller spreads things out
    lf : float
      initial scaling of rank
    tf : int
      total steps
    do_log : bool
    log_n : int
      how many times to output log over the course of time steps
  
  Returns
  -------
  rms : numpy.ndarray, shape (n_save + 1, M, 3)
    location of pseudo atoms, over (evenly spaced intervals of) iterations. 
    0th and 1st iteration are right after each other with nothing in between
  rs : numpy.ndarray, shape (tf+1,3)
    sampled points
  ts_save : numpy.ndarray, shape (n_save + 1,3), dtype np.int32
    time points corresponding to rms
  """
  M = rm0.shape[0]
  d = rm0.shape[1]
  if l0 is None: 0.067*M
  if tf is None: tf=200*M
  rms = np.empty((n_save+1,M,d))
  rms[0] = rm0
  rs = np.empty((tf+1,d))
  ts_save = np.empty((n_save+1)).astype(np.int32)

  r_idxs = np.random.choice(map_3d_idx,p=map_3d_flat,size=tf) # precompute
  rm = rm0
  save_idx=0

  for t in range(tf):
    if do_log and (t % (tf//log_n) == 0): print(t)
    r_idx = r_idxs[t]
    r = xyz[r_idx]
    rs[t]=r
    dist2 = ((r - rm)**2).sum(1) # usually sum over xyz. just need relative rank in one time step, so can use dist2 vs dist
    order = dist2.argsort()
    rank = order.argsort().reshape(-1,1)
    l = l0*(lf/l0)**(t/tf)
    e = e0*(ef/e0)**(t/tf)
    rm = rm + e*np.exp(-rank/l)*(r-rm)
    if do_log and ((t % (tf//n_save) == 0) or (t == tf-1)): 
      rms[save_idx] = rm
      ts_save[save_idx] = t
      save_idx+=1
  return rms,rs,ts_save

def trn_wrapper(map_3d,
                threshold,
                M=1000,
                random_seed=None,
                n_save=10,
                e0=0.3,
                ef=0.05,
                l0_factor=0.005,
                lf=0.5,
                tf_factor=8, # typically M*8
                do_log=True,
                log_n=10
                ):
  map_th = map_3d.copy()
  map_th[map_th < threshold] = 0
  rm0,map_3d_flat,map_3d_idx,xyz,coords_1d = trn_rm0(map_th,M=M,random_seed=random_seed)
  rms,rs,ts_save = trn_iterate(rm0,map_3d_flat,map_3d_idx,xyz,n_save=n_save,e0=e0,ef=ef,l0=M*l0_factor,lf=lf,tf=M*tf_factor,do_log=do_log,log_n=log_n)
  map_3d_th_norm = map_3d_flat.reshape(map_3d.shape)
  return rm0,map_3d_th_norm,map_3d_idx,xyz,coords_1d,rms,rs,ts_save
