import numpy as np

def coords_n_by_d(coords_1d=None,N=None,d=3):
  if N is None: 
    assert coords_1d is not None
  elif coords_1d is None:
    assert N is not None
    coords_1d = np.arange(-N//2,N//2)

  if d==2:
    X = np.meshgrid(coords_1d,coords_1d)
  elif d==3:
    X = np.meshgrid(coords_1d,coords_1d,coords_1d)
  coords = np.zeros((X[0].size,d))
  for di in range(d):
    coords[:,di] = X[di].flatten()
  # make compatible with flatten
  if d == 3: coords[:,[0,1]] = coords[:,[1,0]]
  elif d == 2: coords[:,[0,1]] = coords[:,[1,0]]

  return(coords)

def EA_to_R3(phi, theta, psi=None):
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

def get_random_quat(num_pts,method = 'sphere'):
    """
    Get num_pts of unit quaternions with a uniform random distribution.
    :param num_pts: The number of quaternions to return
    : param method: 
      hemisphere: uniform on the 4 hemisphere, with x in [0,1], y,z in [-1,1]
      sphere: uniform on the sphere, with x,y,z in [-1,1]
    :return: Quaternion list of shape [number of quaternion, 4]
    """
    u = np.random.rand(3, num_pts)
    u1, u2, u3 = [u[x] for x in range(3)]

    quat = np.zeros((4, num_pts))
    if method == 'hemisphere':
      angle = np.pi / 2
    elif method == 'sphere':
      angle = 2 * np.pi
    else:
      assert False, 'use hemisphere or sphere'

    quat[0] = np.sqrt(1 - u1) * np.sin(np.pi * u2 / 2)
    quat[1] = np.sqrt(1 - u1) * np.cos(np.pi * u2 / 2)
    quat[2] = np.sqrt(u1) * np.sin(np.pi * u3 / 2)
    quat[3] = np.sqrt(u1) * np.cos(np.pi * u3 / 2)

    return np.transpose(quat)

def quaternion_to_R(q):
  a,b,c,d = q[0], q[1], q[2], q[3]
  R = np.array([
                [a**2+b**2-c**2-d**2 , 2*b*c-2*a*d , 2*b*d+2*a*c],
                [2*b*c+2*a*d , a**2-b**2+c**2-d**2 , 2*c*d-2*a*b],
                [2*b*d-2*a*c , 2*c*d+2*a*b , a**2-b**2-c**2+d**2]
                ])
  return(R)