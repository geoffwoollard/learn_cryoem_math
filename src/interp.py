import numpy as np

def diff(xy0_rot):
  r0 = np.floor(xy0_rot).astype(np.int)
  r1 = r0+1
  fr = xy0_rot - r0
  mfr = 1-fr
  # assert fr.min() >= 0 and mfr.min() >= 0 
  mfx,mfy,mfz = mfr[:,0], mfr[:,1], mfr[:,-1]
  fx,fy,fz = fr[:,0], fr[:,1], fr[:,-1]
  dd000 = mfz * mfy * mfx;
  dd001 = mfz * mfy *  fx;
  dd010 = mfz *  fy * mfx;
  dd011 = mfz *  fy *  fx;
  dd100 =  fz * mfy * mfx;
  dd101 =  fz * mfy *  fx;
  dd110 =  fz *  fy * mfx;
  dd111 =  fz *  fy *  fx;
  dd = np.array([dd000,dd001,dd010,dd011,dd100,dd101,dd110,dd111])
  return(r0,r1,dd)

def interp_vec(F,r0,r1,dd,N):

  r0_idx = r0 + N//2
  r1_idx = r1 + N//2

  under_grid_idx = np.any(r0_idx < 0,axis=1)
  over_grid_idx = np.any(r1_idx >= N,axis=1)
  good_idx = np.logical_and(~under_grid_idx,~over_grid_idx)


  F_3d_interp = np.zeros((N,N,N)).astype(F.dtype)
  count_3d_interp = np.zeros((N,N,N))
  ones = np.ones(N*N)[good_idx]
  F_flat = F.flatten()[good_idx]

  r0_idx_good = r0_idx[good_idx]
  r1_idx_good = r1_idx[good_idx]

  def fill_vec(F_3d_interp,r0_idx_good,r1_idx_good,F_flat_good,dd):
    dd000,dd001,dd010,dd011,dd100,dd101,dd110,dd111 = dd

    F_3d_interp[r0_idx_good[:,0], r0_idx_good[:,1], r0_idx_good[:,-1]] += F_flat_good*dd000 # 000
    F_3d_interp[r1_idx_good[:,0], r0_idx_good[:,1], r0_idx_good[:,-1]] += F_flat_good*dd001 # 001
    F_3d_interp[r0_idx_good[:,0], r1_idx_good[:,1], r0_idx_good[:,-1]] += F_flat_good*dd010 # 010
    F_3d_interp[r1_idx_good[:,0], r1_idx_good[:,1], r0_idx_good[:,-1]] += F_flat_good*dd011 # 011

    F_3d_interp[r0_idx_good[:,0], r0_idx_good[:,1], r1_idx_good[:,-1]] += F_flat_good*dd100 # 100
    F_3d_interp[r1_idx_good[:,0], r0_idx_good[:,1], r1_idx_good[:,-1]] += F_flat_good*dd101 # 101
    F_3d_interp[r0_idx_good[:,0], r1_idx_good[:,1], r1_idx_good[:,-1]] += F_flat_good*dd110 # 110
    F_3d_interp[r1_idx_good[:,0], r1_idx_good[:,1], r1_idx_good[:,-1]] += F_flat_good*dd111 # 111
    return(F_3d_interp)

  F_3d_interp = fill_vec(F_3d_interp,r0_idx_good,r1_idx_good,F_flat,dd[:,good_idx])
  count_3d_interp = fill_vec(count_3d_interp,r0_idx_good,r1_idx_good,ones,dd[:,good_idx])

    
  return(F_3d_interp,count_3d_interp)