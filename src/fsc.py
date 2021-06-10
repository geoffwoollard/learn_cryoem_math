from geometry import cmask_3d

def complex_corr(arr_1,arr_2):
  return (arr_1*arr_2.conj()).sum()

def do_fsc(half_map_A_f,half_map_B_f,shell_rads):
  fscs = np.zeros(len(shell_rads),dtype=half_map_B_f.dtype)
  N = half_map_A_f.shape[0]
  for idx, radius in enumerate(shell_rads):
    shell_mask = cmask_3d(index=(N//2,N//2,N//2),radius=radius,array=np.ones_like(map_r), do_shell=True, shell_thickness=1).astype(np.bool)
    shell_A = half_map_A_f[shell_mask]
    shell_B = half_map_B_f[shell_mask]
    corr = complex_corr(shell_A,shell_B)
    norm = np.sqrt(complex_corr(shell_A,shell_A)*complex_corr(shell_B,shell_B))
    fscs[idx] = corr / norm
  return fscs