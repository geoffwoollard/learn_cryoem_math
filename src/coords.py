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
  if d == 3: coords[:,[0,1]] = coords[:,[1,0]]

  return(coords)

def EA_to_R3 (phi, theta, psi=None):
    """
    Makes a rotation matrix from Z-Y-Z Euler angles.
    maps image coordinates (x,y,0) view coordinates
    See Z_1 Y_2 Z_3 entry in the table "Proper Euler angles" at https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    http://www.gregslabaugh.net/publications/euler.pdf
    """
    R_z  = np.array([[ np.cos(phi), -np.sin(phi),  0],
                    [ np.sin(phi),  np.cos(phi),  0],
                    [          0,           0,  1]])
    R_y  = np.array([[ np.cos(theta),  0,  np.sin(theta)],
                    [0,              1,             0],
                    [-np.sin(theta),  0,  np.cos(theta)]])
    R = np.dot(R_z, R_y)
    if psi is not None and psi != 0:
        R_in = np.array([[ np.cos(psi), -np.sin(psi),  0],
                        [ np.sin(psi),  np.cos(psi),  0],
                        [          0,           0,  1]])
    
        R = np.dot(R, R_in);

    return R

def deg_to_rad(deg): return(deg*np.pi/180)