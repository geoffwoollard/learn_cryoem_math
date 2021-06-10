import numpy as np
from scipy.ndimage import map_coordinates
import coords

def rotate_map_3d(map_3d, rot, order=1, xyz=None):
  """
  rotate 3d voxelized maps
  this involves several steps
  1. coming up with xyz points that describe the voxel. 
    we have a 3 vector for each voxel. this is made by `coords.coords_n_by_d(N=N,d=2)` and matches the way a 3d array is reshaped. 
    it assumes a 3d array in indexed such that map_3d[x,y,z] is the value at voxel [x,y,z]
  2. rotating the points. the rotation matrix `rot` is 3x3 and `rot.dot(xyz[0]) ` rotates the 3 vector. 
    When we do a list of points with shape (n_points,3) we have to match the shapes, so `xyz_rot = (rot.dot(xyz.T)).T `. 
    Don't forget to  centre coordinates before rotating, otherwise if you use non centred points the origin would be at the top corner of the array.
  3. interpolating the map at the new points
    `scipy.ndimage.map_coordinates(input, coordinates)` returns the value of input (interpolated) at coordinates, where coordinates.shape = (3,n_points). 
    Note that this is the transpose of what `coords.coords_n_by_d` returns. 
    Note that the `coordinates` are relative to the indexing where input[0,0,0] is the top right (not the middle). 
    Stricly speaking this returns the value of the non rotated input at the rotated coordinates, and this rotates the frame and keep the input fixed. 
    This is the same as doing the inverse rotation of the object and keeping the frame fixed. (that's why we use rot.T)
  param:
    map_3d : numpy.ndarray, shape (N,N,N)
    rot : numpy.ndarray, shape (3,3)
  return
    map_3d_rot : numpy.ndarray, shape (N,N,N). rotated map_3d
  """
  assert np.unique(map_3d.shape).size == 1, 'map must be cube, not non-cubic rectangular parallelepiped'
  N = map_3d.shape[0]
  if xyz is None:
    xyz = coords.coords_n_by_d(N=N,d=3) # xyz points cooresponding to the voxel coordinates
  xyz_rot = (rot.T.dot(xyz.T) + N//2)
  map_3d_rot = map_coordinates(map_3d,xyz_rot,order=order).reshape(N,N,N) # reshaped coresponding to xyz. TODO if use mask can do custom reshape
  return map_3d_rot

def cmask_3d(index,radius,array,do_shell=False,shell_thickness=1):
  '''
  make a binary circular mask, or ring (variable thickness). 
  '''
  a,b,c = index
  nx0,nx1,nx2 = array.shape
  x0,x1,x2 = np.ogrid[-a:nx0-a,-b:nx1-b,-c:nx1-c]
  r2 = x0*x0+x1*x1+x2*x2
  mask = r2 <= radius*radius
  if do_shell:
    mask_outer = mask
    mask_inner = r2 <= (radius-shell_thickness)*(radius-shell_thickness)
    mask = np.logical_xor(mask_outer,mask_inner)
  return(mask)