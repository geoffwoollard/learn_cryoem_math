import numpy as np
import numba


def ctf_freqs(N, , psize=1.0, d=2):
  """
  :param shape: Shape tuple.
  :param d: Frequency spacing in inverse Å (1 / pixel size).
  :param full: When false, return only unique Fourier half-space for real data. 
  """

  if d == 1:
    freq_pix_1d = np.arange(0,0.5,1/N)
    return(freq_pix_1d*psize)
  elif d == 2:
    freq_pix_1d = np.arange(-0.5,0.5,1/N)
    x,y = np.meshgrid(freq_1d,freq_1d)
    rho = np.sqrt(x**2+y**2)
    angles_rad = np.arctan2(y, x)
    freq_A_2d = rho * psize
    return(freq_A_2d,angles_rad)

@numba.jit(cache=True, nopython=True, nogil=True)
def eval_ctf(s, a, def1, def2, angast=0, phase=0, kv=300, ac=0.1, cs=2.0, bf=0, lp=0):
    """
    # https://github.com/asarnow/pyem/blob/master/pyem/ctf.py
    :param s: Precomputed frequency grid for CTF evaluation.
    :param a: Precomputed frequency grid angles.
    :param def1: 1st prinicipal underfocus distance (Å).
    :param def2: 2nd principal underfocus distance (Å).
    :param angast: Angle of astigmatism (deg) from x-axis to azimuth.
    :param phase: Phase shift (deg).
    :param kv:  Microscope acceleration potential (kV).
    :param ac:  Amplitude contrast in [0, 1.0].
    :param cs:  Spherical aberration (mm).
    :param bf:  B-factor, divided by 4 in exponential, lowpass positive.
    :param lp:  Hard low-pass filter (Å), should usually be Nyquist.
    """
    angast = np.deg2rad(angast)
    kv = kv * 1e3
    cs = cs * 1e7
    lamb = 12.2643247 / np.sqrt(kv * (1. + kv * 0.978466e-6))
    def_avg = -(def1 + def2) * 0.5
    def_dev = -(def1 - def2) * 0.5
    k1 = np.pi / 2. * 2 * lamb
    k2 = np.pi / 2. * cs * lamb**3
    k3 = np.sqrt(1 - ac**2)
    k4 = bf / 4.  # B-factor, follows RELION convention.
    k5 = np.deg2rad(phase)  # Phase shift.
    if lp != 0:  # Hard low- or high-pass.
        s *= s <= (1. / lp)
    s_2 = s**2
    s_4 = s_2**2
    dZ = def_avg + def_dev * (np.cos(2 * (a - angast)))
    gamma = (k1 * dZ * s_2) + (k2 * s_4) - k5
    ctf = -(k3 * np.sin(gamma) - ac*np.cos(gamma))
    if bf != 0:  # Enforce envelope.
        ctf *= np.exp(-k4 * s_2)
    return ctf