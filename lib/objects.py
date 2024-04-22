import numpy as np
import lib.misc as m
import astropy.units as u
from astropy.coordinates import Distance
from astropy.constants import c, m_p


class block:
    """
    Block class defined having the dimension of the 8x8x8 original block structure. 
    As attribute it has: 'vel', 'dens', 'pres', 'temp', 'energy'.
    Each attribute is computed as the mean value of the center 2x2x2 cells.

    
    For each block the following are computed:

    'n_0', 'distance', 'obs_distance', 'opacity', 'z', 'beta', 'gamma' (lorentz factor)

    
    It has also the following properties:

    'E_i' (Initial Energy due to second order fermi acceleration), 'E_f' (Final energy due to emissions), 'n_e', 'blob', 'ray_line'
    """

    def __init__(self, coords, blocksize, temp, dens, pres, energy, vel, ObsCoords, mu_0=1):
        """
        Parameters:
        --
        coords :class:`~tuple`: tridimensional tuple containing the blocks coordinates

        blocksize :class:`~astropy.unit.Quantity`: block size (GridUnit) 
        
        temp :class:`~astropy.unit.Quantity`: block temperature (K)
        
        dens :class:`~astropy.unit.Quantity`: block density (MassUnit DistanceUnit-3)
        
        pres :class:`~astropy.unit.Quantity`: block pressure (MassUnit DistanceUnit-1 TimeUnit-2)
        
        energy :class:`~astropy.unit.Quantity`: block energy (erg)
        
        vel :class:`~tuple(astropy.unit.Quantity)`: tuple of block velocities on x, y, z 
        
        ObsCoords :class:`~np.array`: array of observer coordinates in cartesian reference frame.
        
        mu_0 :class:`~float`: average atomic mass of swept-up material in units of m_p
        """
        self.x = (coords[0].to(u.cm)).value
        self.y = (coords[1].to(u.cm)).value
        self.z = (coords[2].to(u.cm)).value 
        self.center = np.array([self.x, self.y, self.z])
        self.obs = ObsCoords
        self.radius = (0.5*blocksize.to('cm'))
        self.distance = np.sqrt(coords[0]**2+coords[1]**2+coords[2]**2)
        self.obs_distance = self.BlockObserverDistance()
        self.redshift = self.RedshiftCompiler()

        self.velx = vel[0].to(u.Unit('km s-1'))
        self.vely = vel[1].to(u.Unit('km s-1'))
        self.velz = vel[2].to(u.Unit('km s-1'))
        self.vel = np.sqrt(self.velx**2+self.vely**2+self.velz**2)

        self.temp = temp.to(u.Unit('K'))
        self.dens = dens.to(u.Unit('g cm-3'))
        self.n_0 = (dens/(mu_0*m_p)).to(u.Unit('cm-3')) 
        self.pres = pres.to(u.Unit(' g cm-1 s-2'))
        self.energy = (energy*u.erg)
        self.k = m.opacity(self.temp, self.n_0).k

        self.set_doppler_shift()

        self.CMB = False
        self.EC = False

    def BlockObserverDistance(self):
        """
        Returns the distance of the block from the observer.
        """
        obs_x = self.obs[0]
        obs_y = self.obs[1]
        obs_z = self.obs[2]
        return np.sqrt((self.x-obs_x)**2+(self.y-obs_y)**2+(self.z-obs_z)**2)*u.cm
    
    def ViewingAngle(self):
        """
        Returns the viewing angle wrt to the line of view.
        """
        vel = np.array([self.velx.value, self.vely.value, self.velz.value])

        magn_vel = np.linalg.norm(vel)
        magn_obs = np.linalg.norm(self.obs)

        return np.arccos(np.dot(vel,self.obs)/(magn_vel*magn_obs))
    
    def RedshiftCompiler(self):
        """
        Returns the redshift of the block.
        """
        return Distance(self.obs_distance, unit=self.obs_distance.unit).z
    
    def set_doppler_shift(self):
        """
        Sets the block velocity/c ratio 'beta', doppler shift 'delta_D' annd the lorentz factor 'gamma'
        """
        vel = np.array([self.velx.value, self.vely.value, self.velz.value])
        vel = np.linalg.norm(np.dot(self.obs, vel)/(np.linalg.norm(self.obs)**2)*self.obs)*self.velx.unit
        #vel = 298000*u.Unit('km s-1')
        beta = (vel/c).to(u.Unit(''))
        self.beta = beta.value
        self.delta_D = np.sqrt((1+beta)/(1-beta.value))
        self.gamma = 1/np.sqrt(1-beta.value**2)



    def initial_energy(self, E_i):
        """
        Sets the initial block energy as computed by second order fermi acceleration.

        Parameters:
        --
        E_i: :class:`~float` initial energy in m_e*c^2 units.
        """
        self.E_i = E_i
    
    def final_energy(self, E_f):
        """
        Sets the final block energy computed by taking into account the energy loss due to emissions processes.
        
        Parameters:
        --
        E_f: :class:`~float` final energy in m_e*c^2 units.
        """
        self.E_f = E_f

    def electronic_distribution(self, k_e, n_e, n_e_value):
        """
        Sets the electronic distribution as computed using energy.K_eNormalizer() and agnpy.PowerLaw()

        Paramaters:
        --
        k_e :class:`~astropy.unit.Quantity`: Power law distribution normalization coefficient (m-3)
        
        n_e :class:`~agnpy.PowerLaw`: block electrons distribution function
        
        n_e_value :class:`~astropy.unit.Quantity`: block electron distribution function computed in (E_f+E_i)/2
        """
        self.K_e = k_e
        self.n_e = n_e
        self.n_e_value = n_e_value

    def rayline(self, raypath, target):
        """
        Sets the raypath of the block. 
        
        Parameters:
        --
        raypath :class:`~astropy.unit.Quantity`: thickness of intermediate medium between the source and the observer 
        target :class:`~str`: One end of the line of view, the other being the center of the blob. May be: 'source', 'observer', 'CMB'
        """

        if target == 'source':
            self.src_raypath = raypath
        if target == 'observer':
            self.obs_raypath = raypath
        if target == 'CMB':
            self.cmb_raypath = raypath
    
    def set_blob(self, blob):
        """
        Sets the blob property of the block.

        Parameters:
        --
        blob :class:`~agnpy.emitting_region.Blob`: blob defined by the block characteristics.
        """
        self.blob = blob
    def set_CMB(self, stat=True):
        self.CMB = stat
    def set_EC(self, stat = True):
        self.EC = stat

    def opacity_ssa(self, p, nu, B):
        self.k_ssa = m.ssa_opacity(p, nu, self.K_e, B)
        return self.k_ssa

class target:
    """
    The target to be provided to ray_tracing class. Defined by name and coords.
    """
    def __init__(self, name, coords):
        """
        Parameters:
        --
        name :class:`~str`: name of the target. May be: 'source', 'observer', 'CMB'
        coords :class:`~np.array`: numpy array of the coordinates in a cartesian reference frame.
        """
        self.name = name
        self.coords = coords
