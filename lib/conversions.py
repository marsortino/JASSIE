import numpy as np
import astropy.units as u
from astropy.constants import G, c

def B_to_gaussian_units(B):
    """
    Converts the magnetic field in gaussian units.
    """
    
    B = B.to(u.Unit('G'))
    B = np.sqrt((B.value)**2*u.Unit('cm-1 g s-2'))
    return B

def R_g_units(r, M):
    R_g = G*M/c **2
    return (r / R_g).to('')

def spherical_to_cartesian(rho, theta, phi):
    """
    Converts from spherical to cartesian.
    """
    x = (rho*np.sin(phi)*np.cos(theta)).to('cm')
    y = (rho*np.sin(phi)*np.sin(theta)).to('cm')
    z = (rho*np.cos(phi)).to('cm')
    return np.array([x.value,y.value,z.value])

def angle_between_vectors(vec1, vec2):
    """
    Returns the angle between two vectors.
    vec1 :class: np.array
    vec2 :class: np.array
    """

    magn_vec1 = np.linalg.norm(vec1)
    magn_vec2 = np.linalg.norm(vec2)
    
    return np.arccos(np.dot(vec1,vec2)/(magn_vec1*magn_vec2))