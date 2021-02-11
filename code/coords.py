import numpy as np

def coords_n_by_d(coords_1d=None,N=None,d=3):
  if N is None: 
    assert coords_1d is not None
  elif coords_1d is None:
    assert N is not None
    coords_1d = np.arange(-N//2,N//2)

  # if d==3:
  #   x,y,z = np.meshgrid(coords_1d,coords_1d,coords_1d)
  #   coords = np.zeros((x.size,3))
  #   coords[:,0] = x.flatten()
  #   coords[:,1] = y.flatten()
  #   coords[:,2] = z.flatten()
  if d==2:
    X = np.meshgrid(coords_1d,coords_1d)
  elif d==3:
    X = np.meshgrid(coords_1d,coords_1d,coords_1d)
  coords = np.zeros((X[0].size,d))
  for di in range(d):
    coords[:,di] = X[di].flatten()

  return(coords)