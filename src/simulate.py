from . import coords, interp, fourier, twod, transfer
import mrcfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
import numba

def simulate(map_r,psize,n_particles,snr,N_crop,
		    df_min=15000,
		    df_max=20000,
		    df_diff_min=100,
		    df_diff_max=500,
		    df_ang_min=0,
		    df_ang_max=360,
		    kv=300,
		    cs=2.0,
		    ac=0.1,
		    phase=0,
		    bf=0,
		    do_log=True,
			random_seed=0
			):
	assert np.unique(map_r.shape).size == 1
	N = map_r.shape[0]
	assert N%2 == 0, 'even pixel length'
	
	map_f = fourier.do_fft(map_r)

	N = map_f.shape[0]
	xyz = coords.coords_n_by_d(np.arange(-N//2,N//2),d=3)
	idx_z0 = xyz[:,-1] == 0 
	xy0 = xyz[idx_z0]


	np.random.seed(random_seed)
	qs = coords.get_random_quat(n_particles)
	Rs = coords.quaternion_to_R(qs.T)

	CTFs = transfer.random_ctfs(N,
							    psize,
							    df_min=df_min,
							    df_max=df_max,
							    df_diff_min=df_diff_min,
							    df_diff_max=df_diff_max,
							    df_ang_min=df_ang_min,
							    df_ang_max=df_ang_max,
							    kv=kv,
							    cs=cs,
							    ac=ac,
							    phase=phase,
							    bf=bf,
							    do_log=do_log
							    )


    proj_f = np.zeros((n_particles,N,N),dtype=np.complex64)
	snr = 1
	for idx in range(n_particles):
	  if idx % max(1,(n_particles//10)) == 0: print(idx)
	  R = Rs[:,:,idx]
	  xy0_rot = R.dot(xy0.T).T
	  proj_f[idx] = (map_coordinates(map_f.real, xy0_rot.T + N//2,order=1).astype(np.complex64) + 1j*map_coordinates(map_f.imag, xy0_rot.T + N//2,order=1).astype(np.complex64)).reshape(N,N) # important to keep order=1 for speed. linear is good enough
	  proj_f[idx] *= CTFs[idx]

	i,f = N//2-N_crop//2, N//2+N_crop//2
	proj_r = np.zeros((n_particles,N_crop,N_crop))
	for idx in range(n_particles):
	  proj_r[idx] = twod.do_ifft(proj_f[idx,i:f,i:f]).real
	psize_crop = psize_original*N/N_crop

	signal = np.std(proj_r)
	noise = signal/snr
	proj_r_noise = np.random.normal(loc=proj_r,scale=noise)

	meta_data_df = pd.DataFrame({'df1_A':df1s,'df2_A':df2s,'df_ang_deg':df_ang_deg,'kev':kv,'ac':ac,'cs_mm':cs,
              'rotation_quaternion':[np.array2string(q) for q in qs]
            })

	return(proj_r,proj_r_noise,meta_data_df)