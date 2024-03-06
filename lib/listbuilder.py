import numpy as np
import astropy.units as u
from agnpy.emission_regions import Blob
from lib.objects import block


class builder:
    """
    This class produces a list of objects.block 
    """
    def __init__(self, config):
        """
        Parameters:
        --
        config: :class:`~dictionary` dictionary containing the values needed;  
        """
        self.DistanceUnit = config['DistanceUnit']
        self.TimeUnit = config['TimeUnit']
        self.MassUnit = config['MassUnit']
        self.GridUnit = config['GridUnit']
        self.ObsCoords = config['ObsCoords']
        self.mu_0 = config['mu_0']
        self.B = config['B']

    def build(self, i, data):
        """
        Returns the objects.block, provided the index of the block to be built and the data.

        Parameters:
        --
        i :class:`~int` index of the block to be built;
        data: :class:`~misc.datamap` file data obtained by misc.datamap
        """
        Coords = (data.coordSet[i,0]*self.GridUnit, data.coordSet[i,1]*self.GridUnit, data.coordSet[i,2]*self.GridUnit)
        BlockSize = np.mean(data.BlockSize[i,:])*self.GridUnit
        temp = np.mean(data.TempSet[i, 2:4,2:4,2:4])*u.K
        dens = np.mean(data.DensSet[i, 2:4,2:4,2:4])*(self.MassUnit/self.DistanceUnit**3)
        pres = np.mean(data.PresSet[i, 2:4,2:4,2:4])*self.MassUnit/(self.DistanceUnit*self.TimeUnit**2)
        energy = np.mean(data.EnerSet[i, 2:4,2:4,2:4])*u.erg
        velx = np.mean(data.velxSet[i,2:4,2:4,2:4])*self.DistanceUnit/self.TimeUnit
        vely = np.mean(data.velySet[i,2:4,2:4,2:4])*self.DistanceUnit/self.TimeUnit
        velz = np.mean(data.velzSet[i,2:4,2:4,2:4])*self.DistanceUnit/self.TimeUnit
        vel = (velx, vely, velz)

        data = block(
            Coords,
            BlockSize,
            temp,
            dens,
            pres,
            energy,
            vel,
            self.ObsCoords,
            self.mu_0
        )   

        return data

    def blockbuilder(self, index, data):
        """
        Returns a list of blocks.

        Parameters:
        --
        i :class:`~int` index of the block to be built;
        data: :class:`~misc.datamap` file data obtained by misc.datamap
        """
        blockslist = []
        for i in index:
            blockslist.append(self.build(i, data))
        return blockslist

    def blobbuilder(self, blocklist):
        """
        Given a list of blocks, defines the .blob attribute for each block.
        
        Parameters:
        --
        blocklist :class:`~list` list of objects.block;
        """
        # fopen = open('blobs_properties.txt', 'x')
        # fopen.write('radius       z       delta_D       gamma       B       n_e')
        for block in blocklist:
        #    fopen.write(str(block.radius) + "       " + str(block.z) + "       " + str(block.delta_D) + "       " + str(block.gamma) + "       " + str(self.B) + "       " + str(block.n_e) + "\n")
            blob = Blob(block.radius, block.redshift, block.delta_D, block.gamma, self.B, n_e = block.n_e)
            block.set_blob(blob)
        #fopen.close()
        
    def pointsbuilder(self, index, coordSet):
        """
        Returns an array of points containing the coordinate of each block multiplied by the gridunit. 
        """
        pointlist = np.array([])
        status = 0
        for i in index:
            if status == 0:
                Coords = [coordSet[i,0]*self.GridUnit.value, coordSet[i,1]*self.GridUnit.value, coordSet[i,2]*self.GridUnit.value]
                pointlist = np.array([Coords])
                status += 1
            else:
                Coords = [coordSet[i,0]*self.GridUnit.value, coordSet[i,1]*self.GridUnit.value, coordSet[i,2]*self.GridUnit.value]
                pointlist = np.vstack([pointlist, Coords])
        return pointlist