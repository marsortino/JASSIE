# Misc Utilities 
import h5py
import numpy as np
import sys
import astropy.units as u
from astropy.constants import e, m_e, c, h, k_B, sigma_T
from scipy.special import gamma, factorial 
import lib.conversions as cs
from agnpy.utils.math import phi_to_integrate
from agnpy.utils.conversion import to_R_g_units

import os
import datetime

from decimal import getcontext, Decimal

import lib.quickhull as qh

getcontext().prec = 50
# Order of classes:
# - datamap
# - BlockBuilder
# - BlobBuilder
# - Opacity

class datamap:
    """
    Loads the data in the input file provided in the proper form.
    """

    def __init__(self, config):
        """
        Parameters:
        --
        config :class:`~dictionary`: dictionary containing the main settings.
        """
        self.hdf5data(config)

    def hdf5data(self, config):
        """
        Open the data file and get the main variables
        """
        try:
            fname = h5py.File(config['file_path'], 'r')
        except FileNotFoundError:
            sys.stderr.write("Error: "+ config['file_path']+ ' not found.\n')
            sys.exit(1)
        self.units = (config['DistanceUnit'], config['TimeUnit'], config['MassUnit'])
        self.gridunits = config['GridUnit']
        self.TempSet = np.array(fname.get('temp'), dtype = float)
        self.DensSet = np.array(fname.get('dens'), dtype = float)
        self.PresSet = np.array(fname.get('pres'), dtype = float)
        self.EnerSet = np.array(fname.get('ener'), dtype = float)
        self.velxSet = np.array(fname.get('velx'), dtype = float)
        self.velySet = np.array(fname.get('vely'), dtype = float)
        self.velzSet = np.array(fname.get('velz'), dtype = float)
        self.BlockSize = np.array(fname.get('block size'), dtype = float)
        self.coordSet = np.array(fname.get('coordinates'), dtype = float)
        self.NodeLevel= np.array(fname.get('node type'), dtype=int)

        # Looking for dt min
        self.dtmin = self.timestep(fname['real scalars'])*self.units[1]
        fname.close()
        self.DistanceSet = self.DistanceFromSource(config['GridUnit'])
        self.ObserverDistanceSet  = self.DistanceFromObserver(config['GridUnit'], config['ObserverDistance'])

    def DistanceFromSource(self, GridUnit):
        """
        Returns the distance from the source.
                distance = sqrt(x^2+y^2+z^2)
        
        Parameters:
        --
        GridUnit :class:`~astropy.units.Quantity` Units of the grid in the original simulation box.   
        """
        return np.sqrt(self.coordSet[:,0]**2+self.coordSet[:,1]**2+self.coordSet[:,2]**2)*GridUnit#*self.units[0]
        
    def DistanceFromObserver(self, GridUnit, ObserverDistance = 1e27 * u.cm):
        """
        Returns the observer distance from each blob supposing a edge-on jet.
                        osberver_blob_distance = sqrt(observerdistance**2+z**2)       
        Parameters:
        --
        GridUnit :class:`~astropy.units.Quantity` Units of the grid in the original simulation box.   
        ObserverDistance :class:`~astropy.units.Quantity` Distance of the source from the observer.
        """

        ObserverDistance = ObserverDistance.to(u.cm)
        self.z = self.coordSet[:,2]*GridUnit
        self.y = self.coordSet[:,1]*GridUnit
        self.x = self.coordSet[:,0]*GridUnit
        return (np.sqrt(self.x**2+(ObserverDistance-self.y)**2 + self.z**2)).to(u.cm)

    def timestep(self, parameterpointer):
        """
        Looks for and returns dtmin value on 'real runtime parameters' dataset.
        
        Parameters:
        --
        parameterpointer: pointer to the dataset 'real runtime parameters'.
        """
        i = 0
        for line in parameterpointer:
            a = line[0].decode().replace(" ","")
            if a == 'dt':
                break
            else:
                i+=1 
        return float(parameterpointer[i][1])


class tinymap:
        def __init__(self, config):
            self.read_data(config)
        def read_data(self, config):
            try:
                fname = h5py.File(config['file_path'], 'r')
            except FileNotFoundError:
                sys.stderr.write("Error: "+ config['file_path']+ ' not found.\n')
                sys.exit(1)
            self.units = (config['DistanceUnit'], config['TimeUnit'], config['MassUnit'])
            self.TempSet = np.array(fname.get('temp'), dtype = float)
            self.DensSet = np.array(fname.get('dens'), dtype = float)
            self.PresSet = np.array(fname.get('pres'), dtype = float)
            self.EnerSet = np.array(fname.get('ener'), dtype = float)

class opacity:
    """
    A class that gives back Kramers Opacity for a free-free absorption.
    It gives back the proper opacity value. 
    Opacity is computed as :

                  k_v = \sigma_c * n_e / \rho 

        [ "Astrophysical Concepts - M. Harwit (Third Edition) pp278" - (7-72)]

    where \sigma_c is computed following [Finke 2016] equation (10-11):

                \sigma_c(x) = \sigma_t * S_0(x)
    
    where x is the photon energy in m_e c^2 units and S_0(x) is:

                S_0(x) = 3/(8x^2)*[4+ (2x^2*(1+x))/(1+2x)^2) + (x^2-2x-2)/x*ln(1+2x)]

    """
    def __init__(self, temp, n_0):
        """
        Parameters:
        --
        temp : :class:`~astropy.units.Quantity` temperature value of the block;
        n_0 : :class:`~astropy.units.Quantity` number density value of the block;
        """
        nu = 4.965*k_B*temp/h
        a = (h*nu/(m_e*c**2)).to(u.Unit(""))
        self.k = self.ComputeOpacity(n_0, a)

    def ComputeCrossSection(self, x):
        """
        Computes Total Compton Cross Section, calculating S_0(x) from equation (11) [Finke 2016]
        Parameters:
        --
        x : :class:`~astropy.units.Quantity`: photon energy in m_e c^2 units.
        """
        prefactor = 3/(8*x**2)
        factor = (4+ ((2*x**2)*(1+x))/(1+2*x)**2 + (x**2 -2*x -2)/(x)*np.log(1+2*x))
        return sigma_T * prefactor * factor


    def ComputeOpacity(self, n_0, x):
        """
        Opacity is computed as :

            k = \sigma * n_e / \rho 

            i.e.
                k = \sigma *n_0
        [ "Astrophysical Concepts - M. Harwit (Third Edition) pp278" - (7-72)]

        Parameters:
        --
        n_0 : :class:`~astropy.units.Quantity` number density value of the block;
        x : :class:`~astropy.units.Quantity`: photon energy in m_e c^2 units.  

        """
        sigma_c = self.ComputeCrossSection(x)
        return sigma_c * n_0

def mu_from_r_tilde_dec(R_in_tilde, R_out_tilde, r_tilde, size = 100):
    r"""array of cosine angles, spanning from :math:`R_{\mathrm{in}}` to
    :math:`R_{\mathrm{out}}`, viewed from a given height :math:`\tilde{r}`
    above the disk, Eq. 72 and 73 in [Finke2016]_.
    
    From agnpy module. Modified since sometimes a mantissa of more than 15 values is needed..
    N.b. np.linspace does not work with Decimal object type.
    """


    mu_min = 1 / np.sqrt(1 + np.power((R_out_tilde / r_tilde), 2))
    mu_max = 1 / np.sqrt(1 + np.power((R_in_tilde / r_tilde), 2))

    output = np.array([mu_min + i*(mu_max-mu_min)/size for i in range(size)])

    return output


# def compare():
class Compare:
    """
    Reduces the original index after removing the blocks that not compose the sphere or half sphere in cmb_blocks and disc_blocks
    """
    def __init__(self, index_size):
        self.zeros = np.zeros(index_size)
        
    def IndexReduce(self, tiny_index):
        j = 0
        i = 0
        for value in tiny_index:
            while j != value:
                i += 1
                if self.zeros[i] == 0:
                    j += 1

            while self.zeros[i] != 0:
                i += 1
            self.zeros[i] = 1

    def blocks_selection(self, blocklist, case):
        """
        Set the CMB or EC properties of the block object as True. 
        """
        index = np.where(self.zeros == 1)[0]
        if case == 'CMB':
            for value in index:
                blocklist[value].set_CMB()
        if case == 'EC':
            for value in index:
                blocklist[value].set_EC()

    
def sphere(pointlist, blocklist, qhull_depth):
    """
    Computes the subset of blocks which compose the external shells of the distribution of blocks.

    Parameters:
    --
    pointlist :class:`~list`: list of blobs coordinates.
    blocklist :class:`~list[objects.block]`: list of blocks.
    qhull_depth :class:`~float`: value ranging from 0 to 1 indicating the percentage of blocks to consider.
     """
    sec_index = 0
    lenght_pointlist = len(pointlist)*qhull_depth
    reduce_index = Compare(len(pointlist))
    while (sec_index <= lenght_pointlist):
        points_index = qh.convex_hull(pointlist)
        reduce_index.IndexReduce(points_index)
        point_index = np.flip(points_index)
        pointlist = np.delete(pointlist, point_index, 0)
        sec_index += len(points_index)

    reduce_index.blocks_selection(blocklist, 'CMB')


def get_lower(polygon):
    """
    Returns the blocks whose whose z-axis value is below the midpoint of the distribution's z-axis

    Parameters:
    --
    polygon :class:`~list` list of points indicating the vertices of the polygon
    """
    minz = np.min(polygon[:,2])
    maxz = np.max(polygon[:,2])
    if abs(minz)>abs(maxz):
        # If the sphere is on the negative z-axis:
        tmp = minz
        minz = maxz
        maxz = tmp
        differenza = (maxz-minz)
        differenza = differenza*0.01+minz
        lower_curve = np.where(polygon[:,2]>differenza)[0]
        #lower_curve = np.logical_and(polygon[:,2]>differenza, polygon[:,2]<threeshold)
        return lower_curve
    
    differenza = (maxz-minz)
    differenza = differenza*0.01+minz
    lower_curve = np.where(polygon[:,2]<differenza)[0]
    return lower_curve

def lower_sphere(pointlist, blocklist, qhull_depth):
    """
    Computes the subset of blocks which compose the external shells of the distribution of blocks. Only the blocks near the source are considered.

    Parameters:
    --
    pointlist :class:`~list`: list of blobs coordinates.
    blocklist :class:`~list[objects.block]`: list of blocks.
    qhull_depth :class:`~float`: value ranging from 0 to 1 indicating the percentage of blocks to consider.
    """
    sec_index = 0
    lenght_pointlist = len(pointlist)*qhull_depth/2.0 #Since we are taking the lower part of the sphere, the index must have length/2 wrt to whole sphere index lenght
    reduce_index_disc = Compare(len(pointlist))
    while(sec_index <= lenght_pointlist):
        points_vertices = qh.convex_hull_vertices(pointlist)
        points_vertices = get_lower(points_vertices)
        reduce_index_disc.IndexReduce(points_vertices)
        point_vertices = np.flip(points_vertices)
        pointlist = np.delete(pointlist, point_vertices, 0)
        sec_index += len(points_vertices)
    reduce_index_disc.blocks_selection(blocklist, 'EC')


class integrals:
    """
    Computes integrals using the trapezoidal method compatible with the Decimal routine.
    
    provvisory. to be removed.
    
    """
    def __init__(self, R_in, r) -> None:
        self.R_in = R_in
        self.r = r

    def func(self, x):
        x = x**(-2)-1
        sqrt = pow(x, Decimal('0.5'))
        num = 1 - pow(self.R_in/(self.r*sqrt), Decimal('0.5'))
        den = pow(x, Decimal('1.5'))
        return num/den

    def trapz(self, x):
        n = len(x)
        a = x[0]
        b = x[n-1]
        h = (b - a) / Decimal(n)
        integral_sum = Decimal('0.5') * (self.func(a) + self.func(b))

        for i in range(n):
            x_i = a + Decimal(i) * h
            integral_sum += self.func(x_i)

        result = h * integral_sum
        return result


def ssa_coefficient_c(p):
    """
    Gould 1979 eq 32
    """
    prefactor = pow(3, (p+1)/2)/8
    gamma_1 = gamma((p+6)/4)
    gamma_2 = gamma((3*p+2)/12)
    gamma_3 = gamma((3*p+22)/12)
    gamma_4 = gamma((p+8)/4)

    return prefactor*pow(np.pi, 0.5)*gamma_1*gamma_2*gamma_3/gamma_4
def ssa_opacity(p, nu, K_e, B):
    """
    Dermer Menon 7.148

    p coefficiente power law
    nu frequenza
    K_e constante normalizzazione pl
    nu_B gyrofrequency
    c velocità luce
    c(p) coefficiente per opacity

    k_nu = c(p) \frac{c r_e K_e}/nu (nu_B/nu)^((p+2)/2)
    """

    B = cs.B_to_gaussian_units(B)
    coeff = ssa_coefficient_c(p)
    r_e = (e.esu**2/(m_e*c**2)).to(u.cm)
    nu_B = (e.esu*B/(2*np.pi*m_e*c)).to(u.Hz)
    prefactor = (c*r_e*K_e/nu).to(u.Unit('cm-1'))
    exp = ((nu_B/nu).to(u.Unit('')))**((p+2)/2)

    return coeff*prefactor*exp
    

def BlocklistTXT(blocklist, filename):
    """
    Write a txt file containing:
        First Row: the number of the blocks
        Following Rows: x, y, z and radius of each block
    
    Each value is adimensional. 

        x,y,z are in cm, radius is in cm
    
    """

    file_path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(file_path, 'tmp_files')
    if not os.path.exists(file_path):
        os.mkdir(file_path)
    
    file_path = os.path.join(file_path, 'tmp_list_value_' + filename+ '.txt')

    if os.path.exists(file_path):
        block_list = open(file_path, 'w')
    else:
        block_list = open(file_path, 'x')

    block_list.write(str(len(blocklist)))
    for block in blocklist:
        block_list.write('\n'+str(block.x)+' '+str(block.y)+' '+str(block.z)+' '+str(block.radius.value))

    block_list.close()


def write_log(config, time_table, n_blocks, final_flux):
    """
    Writes a final log containing the parameters used for the computation.
    
    Parameters:
    --
    config :class: `~dictionary`: contains all initial settings.
    time_table :class: `~dictionary`: contains total time for each computation. 
    n_blocks :class: `~int`: total number of blocks
    """
    log_dir = 'logs/'

    final_log_name = 'log_'+config['name']+'.txt'
    final_log_name = os.path.join(log_dir, final_log_name)

    key_list = {
        "Launch Settings": ('id_launch', 'num_cores', 'qhull_depth'),
        "Jet Settings": ('n_e', 'p_min', 'p_max', 'gamma_Min', 'B', 'Z', 'mu_0', 'E_released', 'f_delta'),
        "Plot Settings": ('nu_Min', 'nu_Max', 'ObserverDistance', 'Obs_theta', 'Obs_phi'),
        "Unit Settings": ('DistanceUnit', 'TimeUnit', 'MassUnit', 'GridUnit'),
    }
    
    disc_settings = [
        'L_disc',
        'eta',
        'R_g_unit',
        'R_in',
        'R_out',
    ]


    f_log = open(final_log_name, 'x')
    f_log.write('Recap of file ' + str(config['name']) + ':\n')
    
    if config['blocks'] == 'all':
        f_log.write('Blocks selected: all\n')
    elif config['blocks'] == 'both':
        write_blocks_both(f_log, config)
    elif config['blocks'] == 'phys_only':
        write_log_physical_quantities(f_log, config)
    elif config['blocks'] == 'coords_only':
        write_log_physical_quantities(f_log, config)
    f_log.write('for a total number of blocks of ' + str(n_blocks) + '.\n')
    f_log.write('Processes considered:\n')
    if config['id_brem']:
        f_log.write('-- Bremmstralung\n')
    if config['id_syn']:
        f_log.write('-- Synchrotron\n')
        if config['id_ssa']:
            f_log.write('-- -- Self Absorption\n')
        if config['id_ssc']:
            f_log.write('-- -- Self Compton\n')
    if config['id_ec']:
        f_log.write('-- External Compton\n')
        if config['id_cmb']:
            f_log.write('-- -- CMB\n')
        if config['id_disc']:
            f_log.write('-- -- Disc\nDisc Settings:\n')
            
            if 'M_BH' in config:
                f_log.write('M_BH: ' + str(config['M_BH']) + '\n')
            else:
                f_log.write('a: ' + str(config['a']) + '\n')
                f_log.write('P_jet: ' + str(config['P_jet']) + '\n')

            for key in disc_settings:
                f_log.write(str(key)+": " + str(config[key]) + '\n')

    for key in key_list:
        f_log.write(str(key) + '\n')
        for input in key_list[key]:
            f_log.write(str(input)+': '+ str(config[input]) + '\n')

    f_log.write('Timings:\n')
    for key in time_table:
        f_log.write('-- ' +key + ': ' + str(time_table[key]) + '\n')

    i = 0
    f_log.write('Final Flux (in '+ str(final_flux[0].unit)+ '):\n')
    for sed_value in final_flux:
        f_log.write(str(sed_value.value))
        if i == 5:
            f_log.write('\n')
            i = 0
        i+=1

    f_log.write('-----------------------')
    date = datetime.datetime.now()
    f_log.write("Current time: %s/%s/%s, %s:%s:%s" %(date.year, date.month, date.day, date.hour, date.minute, date.second))
    f_log.write('-----------------------')
    f_log.close()

def clean(name):
    """
    Removes some temporary files.
    """
    hidden_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'logs/tmp/')
    hidden_file_name = os.path.join(hidden_dir, 'tmp_' + name)
    try:
        os.remove(hidden_file_name)
        print('Cleaned.')
    except FileNotFoundError:
        print('Cannot find ' +'tmp_'+name+' which is strange.')
        print('Code may have still worked, however be wary of overwriting issues if different setups of same input file have been launched together.')

def write_blocks_both(file_pointer, config):
    """
    Writes on final log files both physical quantities and coordinates used for blocks selecting
    """

    quantity = [
        'x_Min',
        'x_Max',
        'y_Min',
        'y_Max',
        'z_Min',
        'z_Max',
        'pres_Min',
        'pres_Max',
        'dens_Min',
        'dens_Max',
        'energy_Min',
        'energy_Max',
        'T_Min',
        'T_Max',
    ]

    file_pointer.write('Blocks selected both physical quantities and coordinates:\n')
    
    for parameter in quantity:
        file_pointer.write(parameter+': '+str(config[parameter])+'\n')

def write_log_physical_quantities(file_pointer, config):
    """
    Writes on final log files physical quantities used for blocks selecting

    """
    quantity = [
        'pres_Min',
        'pres_Max',
        'dens_Min',
        'dens_Max',
        'energy_Min',
        'energy_Max',
        'T_Min',
        'T_Max',
    ]
    file_pointer.write('Blocks selected based on physical quantities:\n')
    for parameter in quantity:
        file_pointer.write(parameter+': '+str(config[parameter])+'\n')

def write_log_coords(file_pointer, config):
    """
    Writes on final log files coordinates used for blocks selecting
    """

    quantity = [
        'x_Min',
        'x_Max',
        'y_Min',
        'y_Max',
        'z_Min',
        'z_Max',
    ]

    file_pointer.write('Blocks selected based on coordinates:\n')
    for parameter in quantity:
        file_pointer.write(parameter+': '+str(config[parameter])+'\n')

def get_linear_gradient(p_min, p_max, blocklist):
    """
    """
    n_block = len(blocklist)

    distance_min = np.abs(blocklist[0].distance)
    distance_max = distance_min+1*distance_min.unit

    for i in range(1, n_block):
            block_distance = blocklist[i].distance
            if distance_min > block_distance:
                distance_min = block_distance
            elif distance_max < block_distance:
                distance_max = block_distance

    print('min distance from center:', distance_min)
    print('max distance from center:', distance_max)

    p_array = []
    for i in range(0, n_block):
        block_distance = blocklist[i].distance
        p_array.append(p_min + (p_max - p_max)*(block_distance - distance_min)/(distance_max - distance_min))
        

    return p_array