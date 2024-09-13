import sys
import numpy as np
import astropy.units as u
from astropy.constants import m_p


def MainIndex(data, config):
    """
    Returns the index of the blocks that satisfies the conditions included in config:
        i.e.:
                        dens_Min < dens < dens_Max
                        pres_Min < pres < pres_Max
                        temp_Min < temp < temp_Max
                        energy_Min < energy < energy_Max

    Parameters:
    --
    data: :class:`~misc.datamap` file data obtained by misc.datamap
    config: :class:`~dictionary` dictionary containing the values needed;
    """

    if config['blocks'] == 'all':
        Index = np.unique(np.argwhere(data.DistanceSet>0))
        return Index

    elif config['blocks'] == 'phys_only':
        Index = Index_quantities(data, config)

    elif config['blocks'] == 'coords_only':
        Index = Index_coords(data, config)
    
    elif config['blocks'] == 'both':
        Index = Index_both(data, config)

    Index_source =  np.unique(np.argwhere(data.DistanceSet>0))
    Index = np.intersect1d(Index, Index_source)

    if not np.any(Index):
        print('No block satisfies the initial conditions.\n')
        sys.exit(1)

    return Index

def Index_quantities(data, config):
    """
    Returns the index of blocks that satisfy the physical conditions stated in config.txt
    """

    n = len(data.TempSet[:,1,1,1])

    n_0 = np.zeros(n)
    density = np.zeros(n)
    pressure = np.zeros(n)
    energy = np.zeros(n)
    temp = np.zeros(n)
    distance_source = np.zeros(n)

    mu_0 = config['mu_0']
    
    for i in range(n):
        density[i] = (np.mean((data.DensSet[i,2:4,2:4,2:4]*data.units[2]/(data.units[0]**3)).to(u.Unit('g cm-3'))).value)
        pressure[i] = ((np.mean(data.PresSet[i, 2:4,2:4,2:4]*data.units[2]/(data.units[0]*data.units[1]**2)).to(u.Unit('g cm-1 s-2'))).value)
        n_0[i] = (np.mean((((data.DensSet[i,2:4,2:4,2:4])*data.units[2]/(mu_0*m_p*data.units[0]**3)).to(u.Unit('cm-3'))).value)) # Converted in cm^-3
        temp[i] = np.mean(data.TempSet[i,2:4,2:4,2:4])
        energy[i] = (np.mean(data.EnerSet[i, 2:4,2:4,2:4])*u.erg).value
        distance_source[i] = data.DistanceSet[i].value

    IndexDens = np.unique(np.argwhere((density > config['dens_Min'].value) & (density < config['dens_Max'].value)))
    IndexPres = np.unique(np.argwhere((pressure > config['pres_Min'].value)& (pressure < config['pres_Max'].value)))

    IndexTemp = np.unique(np.argwhere((temp> config['T_Min'].value) & (temp < config['T_Max'].value)))
    IndexEnergy = np.unique(np.argwhere((energy > config['energy_Min'].value)& (energy < config['energy_Max'].value)))

    Index1 = np.intersect1d(IndexDens, IndexPres)
    Index2 = np.intersect1d(IndexTemp, IndexEnergy)

    Index = np.intersect1d(Index1, Index2)

    return Index

def Index_coords(data, config):
    """
    Returns the index of blocks that satisfy the spatial conditions stated in config.txt
    """

    IndexCoord1 = np.unique(np.argwhere((data.x >= config['x_min']) & (data.x <= config['x_max'])))
    IndexCoord2 = np.unique(np.argwhere((data.y >= config['y_min']) & (data.x <= config['y_max'])))
    IndexCoord3 = np.unique(np.argwhere((data.z >= config['z_min']) & (data.x <= config['z_max'])))

    Index = np.intersect1d(np.intersect1d(IndexCoord1, IndexCoord2), IndexCoord3)

    return Index

def Index_both(data, config):
    """
    Returns the index of blocks that satisfy both the spatial and physical conditions stated in config.txt
    """

    n = len(data.TempSet[:,1,1,1])
    n_0 = np.zeros(n)
    density = np.zeros(n)
    pressure = np.zeros(n)
    energy = np.zeros(n)
    temp = np.zeros(n)
    distance_source = np.zeros(n)


    mu_0 = config['mu_0']
    
    for i in range(n):
        density[i] = (np.mean((data.DensSet[i,2:4,2:4,2:4]*data.units[2]/(data.units[0]**3)).to(u.Unit('g cm-3'))).value)
        pressure[i] = ((np.mean(data.PresSet[i, 2:4,2:4,2:4]*data.units[2]/(data.units[0]*data.units[1]**2)).to(u.Unit('g cm-1 s-2'))).value)
        n_0[i] = (np.mean((((data.DensSet[i,2:4,2:4,2:4])*data.units[2]/(mu_0*m_p*data.units[0]**3)).to(u.Unit('cm-3'))).value)) # Converted in cm^-3
        temp[i] = np.mean(data.TempSet[i,2:4,2:4,2:4])
        energy[i] = (np.mean(data.EnerSet[i, 2:4,2:4,2:4])*u.erg).value
        distance_source[i] = data.DistanceSet[i].value

    IndexDens = np.unique(np.argwhere((density > config['dens_Min'].value) & (density < config['dens_Max'].value)))
    IndexPres = np.unique(np.argwhere((pressure > config['pres_Min'].value)& (pressure < config['pres_Max'].value)))

    IndexTemp = np.unique(np.argwhere((temp> config['T_Min'].value) & (temp < config['T_Max'].value)))
    IndexEnergy = np.unique(np.argwhere((energy > config['energy_Min'].value)& (energy < config['energy_Max'].value)))

    IndexCoord1 = np.unique(np.argwhere((data.x >= config['x_min']) & (data.x <= config['x_max'])))
    IndexCoord2 = np.unique(np.argwhere((data.y >= config['y_min']) & (data.x <= config['y_max'])))
    IndexCoord3 = np.unique(np.argwhere((data.z >= config['z_min']) & (data.x <= config['z_max'])))

    Index1 = np.intersect1d(IndexDens, IndexPres)
    Index2 = np.intersect1d(IndexTemp, IndexEnergy)
    Index3 = np.intersect1d(np.intersect1d(IndexCoord1, IndexCoord2), IndexCoord3)

    Index = np.intersect1d(np.intersect1d(Index1, Index2), Index3)

    return Index


def SecIndex(Index, n_0, start, stop):
    """
    TBD
    """
    Index1 = np.unique(np.argwhere((start<n_0) & (n_0<stop)))
    if not np.any(Index1):
        print("Found no block with n_0 values between: ", str(start), 'and ', str(stop), '.\n')
        return False
    return np.intersect1d(Index1, Index)

class PlotIndex:
    """
    similar to MainIndex, used in plotting.py.
    Returns the DataSet needed.
    """
    def __init__(self, data, config):
        self.data = data
        self.config = config
        self.plot_index()

    def check_len_zero(self, index, label):
        if len(index) == 0:
            print('The interval choosed for', label, 'is empty. Exiting.')
            sys.exit(1)

    def plot_index(self):

        n = len(self.data.TempSet[:,1,1,1])
        density = np.zeros(n)
        pressure = np.zeros(n)
        energy = np.zeros(n)
        temp = np.zeros(n)

        for i in range(n):
            density[i] = (np.mean((self.data.DensSet[i,2:4,2:4,2:4]*self.data.units[2]/(self.data.units[0]**3)).to(u.Unit('g cm-3'))).value)
            pressure[i] = ((np.mean(self.data.PresSet[i, 2:4,2:4,2:4]*self.data.units[2]/(self.data.units[0]*self.data.units[1]**2)).to(u.Unit('g cm-1 s-2'))).value)
            temp[i] = np.mean(self.data.TempSet[i,2:4,2:4,2:4])
            energy[i] = (np.mean(self.data.EnerSet[i, 2:4,2:4,2:4])*u.erg).value


        Index = np.array(range(0, len(density)))
        if self.config['kind'] == 'histogram':
            for label in self.config['label_list']:
                if label == 'density':
                    index_label = np.unique(np.argwhere((density > self.config['dens_Min'].value) & (density < self.config['dens_Max'].value)))
                elif label == 'pressure':
                    index_label = np.unique(np.argwhere((pressure > self.config['pres_Min'].value)& (pressure < self.config['pres_Max'].value)))
                elif label == 'temperature':
                    index_label = np.unique(np.argwhere((temp> self.config['T_Min'].value) & (temp < self.config['T_Max'].value)))
                elif label == 'energy':
                    index_label = np.unique(np.argwhere((energy > self.config['energy_Min'].value)& (energy < self.config['energy_Max'].value)))
                self.check_len_zero(index_label, label)
                Index = np.intersect1d(Index, index_label)

        elif self.config['kind'] == 'contour':
            label_list = [self.config['x_label'], self.config['y_label']]
            for label in label_list:
                if label == 'density':
                    index_label = np.unique(np.argwhere((density > self.config['dens_Min'].value) & (density < self.config['dens_Max'].value)))
                elif label == 'pressure':
                    index_label = np.unique(np.argwhere((pressure > self.config['pres_Min'].value)& (pressure < self.config['pres_Max'].value)))
                elif label == 'temperature':
                    index_label = np.unique(np.argwhere((temp> self.config['T_Min'].value) & (temp < self.config['T_Max'].value)))
                elif label == 'energy':
                    index_label = np.unique(np.argwhere((energy > self.config['energy_Min'].value)& (energy < self.config['energy_Max'].value)))
                Index = np.intersect1d(Index, index_label)
        
        density =density[Index] 
        pressure = pressure[Index]
        energy = energy[Index]
        temp = temp[Index]

        self.density = density
        self.pressure = pressure
        self.temp = temp
        self.energy = energy
