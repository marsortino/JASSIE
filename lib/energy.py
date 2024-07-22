import numpy as np
import astropy.units as u
from astropy.constants import e, m_e, m_p, c, alpha, sigma_T, mu0, G
from agnpy.targets import SSDisk, RingDustTorus, SphericalShellBLR, CMB
from agnpy.spectra import PowerLaw, BrokenPowerLaw
from agnpy.emission_regions import Blob
import lib.conversions as cs
from lib.misc import mu_from_r_tilde_dec, integrals

from decimal import Decimal, getcontext 

getcontext().prec = 50

class Block_Energy:
    """
    Computes the energy of the block. For each block it computes the maximum electron energy as reached
    via Second Order Fermi Acceleration. 
    The final energy is computed by taking into account the emission processes.
    """
    r_e = (e.esu**2/(m_e*c**2)).to(u.cm)

    def __init__(self, config):
        """
        Parameters:
        --
        config: :class:`~dictionary` dictionary containing the values needed
        external_blocks: :class:`~dictionary` dictionary containing subsets of the external blocks:
                 - positive_sphere: all external blocks on the positive z axis
                 - negative_sphere: all external blocks on the negative z axis
                 - positive_half_sphere: external blocks near the source with positive z axis
                 - negative_half_sphere: external blocks near the source with negative z axis
        """
        self.B = config['B']
        self.Z = config['Z']
        self.mu_0 = config['mu_0']
        if config['n_e'] == 'PowerLaw':
            self.n_e = 'PowerLaw'
            self.p = config['p']
        elif config['n_e'] == 'BrokenPowerLaw':
            self.n_e = 'BrokenPowerLaw'
            self.p1 = config['p1']
            self.p2 = config['p2']
            self.gamma_b = config['gamma_B']
        self.E_released = config['E_released']
        self.f_delta = config['f_delta']
        self.gamma_min = config['gamma_Min']

        self.id_brem = config['id_brem']
        self.id_syn = config['id_syn']

        self.sphere = False


    def CoolingTime(self, block):
        """
        Computes the cooling time. Given by:
                t = E/(-dE/dt)
        
        where (-dE/dt) is given by Blumenthal & Gould (1970) eq. [3.53]:
            -(dE/dt) = 4\alpha r_0^2 * c * [\sum_z n_z Z(Z+1)(ln (2E_i)-1/3)]E_i

        this means that:
                t = {4\alpha r_0^2 * c * [\sum_z n_z Z(Z+1)(ln (2E_i)-1/3)]}^-1
        """

        t = 4*alpha*self.r_e**2*c*block.n_0*self.Z*(self.Z+1)*(np.log(2*block.E_i)-1/3)
        return (1/t).to(u.s)
    
    def SynchroTime(self, block):
        """
        Computes the Cooling Time: 

            t_Syn = m_e*c/sigma_t * (B^2 sin^2 alpha/mu_0)^-1 * gamma_0^-1

        where gamma_0^-1 is from the electron initial energy E_0 = gamma_0 m_e c^2
        """

        t = ((m_e*c/sigma_T*(self.B**2/mu0)**(-1))/block.E_i).to(u.s)
        return t
    
    
    def ExtComptonCoolingTime(self, target, block):
        """
        Returns the external compton cooling time and density radiation field.
        Computes the External compton cooling time for both CMB and SSDisk cases.

            t_EC = 3/4 m_e*c/sigma_T * U_rad^{-1} gammma_0^-1

        where U_rad is the energy density of the photons before being scattered.
        In the disc case U_rad is computed as:
        
            U_rad = 3/(8 \pi c)* (GMm_dot)/r^3 \int_mu_min^mu_max (1-sqrt(R_in/(r*sqrt(mu^-1-1)))) (mu^-2-1)^-3/2 
        
        """
        r = block.distance
        if isinstance(target, CMB):
            U_rad = target.u_0*np.exp(-block.k*block.obs_raypath).to(u.Unit(''))

        if isinstance(target, SSDisk):
            M_BH = target.M_BH
            L_disk = target.L_disk
            eta = target.eta
            R_in = target.R_in
            R_out = target.R_out


            # SUBJECT TO BE CHANGED. MPMATH, TO REMOVE ALSO BLOCKS TOO FAR AWAY.
            r_tilde = cs.R_g_units(r, M_BH).value
            R_in_tilde = cs.R_g_units(R_in, M_BH).value
            R_out_tilde = cs.R_g_units(R_out, M_BH).value
            if np.sqrt(1+np.power(R_in_tilde/r_tilde,2)) == 1:
                block.set_EC(False)
                return None, 0
            m_dot = (L_disk / (eta * np.power(c, 2))).to("g / s")
            # for the disk we do not integrate mu from -1 to 1 but choose the range
            # of zenith angles subtended from a given distance
            mu = mu_from_r_tilde_dec(R_in_tilde, R_out_tilde, r_tilde, 100)
            integrand = (1-np.sqrt(R_in_tilde/(r_tilde*np.sqrt(np.power(mu,-2)-1))))/(np.power(mu**(-2)-1,3/2))
            integral = np.trapz(integrand, mu)

            # class_integral = integrals(R_in_tilde, r_tilde)
            # integral = class_integral.trapz(mu)

            # r_tilde = float(r_tilde)
            # integral = float(integral)
            prefactor = 3*G*M_BH*m_dot/(8*np.pi*c*r**3*np.power(r_tilde,3))
            U_rad = (prefactor*integral).to(u.Unit('erg cm-3'))
            U_rad = U_rad*np.exp((-block.k*block.src_raypath).to(u.Unit(''))) 

        t_EC = ((3/4*(m_e*c/sigma_T*U_rad**(-1)*block.E_i**(-1))).to(u.s))
        return t_EC, U_rad 
    
        
    def E_initial(self, block):
        """
        Computes the initial energy in m_e c**2 units, using Dermer & Menon (2009) eq: 14.90

                        E = 7.7 * 10^20 Z epsilon_B^1/2 (mu_0 n_0)^1/6 f_delta beta_0 (E_released Gamma_0})^1/3 eV
        where f_delta: shell scale of swept up material
        epsilon_B: fraction of total energy of magnetic field 'epsilon_B = B^2/(8\pi e) where e is internal energy
        E_released: apparent isotropy energy release of the explosion in rest mass solar unit                


        Note: while the formula states the part (mu_0 n_0) and the latter is defined in unit of cm-3, all the units of measure are already taken into account in eV.
        Hence n_0 is used.
        """
        block.set_doppler_shift()
        gamma = block.gamma
        dens = block.dens
        beta = block.beta
        n_0 = block.n_0
        B = cs.B_to_gaussian_units(self.B)
        e = (4*gamma**2*dens*c**2).to(u.Unit('g m-1 s-2'))
        epsilon_B = (B**2/(8*np.pi*e))
        E_initial = (7.7*10e20*self.Z*epsilon_B**0.5*(self.mu_0*n_0.value)**(1/6)*self.f_delta*beta*(self.E_released*gamma)**(1/3)).to(u.Unit(''))
        adimensional_unit = (1*u.eV/(m_e * c**2)).to(u.Unit(''))
        block.initial_energy(E_initial.value*adimensional_unit.value)

    def E_final(self, block, dtmin, targetlist = None, cmb = False):
        """
        Computes the final energy in m_e c**2 units of the block after taking into account the energy loss due to emission.
        
        Parameters:
        --
        block :class:`~objects.block`: block to compute.
        dtmin :class:`~astropy.unit.Quantity`: time interval between the current file simulation's time and the following one.  
        targetlist :class:`~list[agnpy.target]`: list of targets that interact with the block via external compton
        cmb :class:`~boolean`: check whether to consider the CMB in the External Compton computation or not
        
        """
        t = []
        u_rad = []

        if self.id_brem:
            t_brem = self.CoolingTime(block)        
            t.append(t_brem)
        else:
            t_brem = 99e99 * u.s

        if self.id_syn:
            t_syn = self.SynchroTime(block)    
            t.append(t_syn)

        else:
            t_syn = 99e99 * u.s

        if targetlist != None:
            if block.EC == True:
                for target in targetlist:
                    ExtCompton_values = self.ExtComptonCoolingTime(target, block)
                    if ExtCompton_values[0] == None:
                        continue
                    t.append(ExtCompton_values[0])
                    u_rad.append(ExtCompton_values[1])

            

        if cmb == True:
            if block.CMB == True:
                ExtCompton_values = self.ExtComptonCoolingTime(CMB(block.redshift), block)
                t.append(ExtCompton_values[0])
                u_rad.append(ExtCompton_values[1])

        t.append(98e99 * u.s)
        t_min = min(t)

        if t_min == 98e99 * u.s:
            E_m = 1
            block.final_energy(E_m)
            return

        elif t_min <= dtmin:
            E_m = self.gamma_min+0.1

        elif t_min == t_brem:
            E_m = (block.E_i*(1-4*alpha*self.r_e**2*c*block.n_0*self.Z*(self.Z+1)*(np.log(2*block.E_i)-1/3)*dtmin)).to(u.Unit('')).value
        
        elif t_min == t_syn:
            E_m = (((-4/3*sigma_T*c*self.B**2/(2*mu0)*block.E_i**2)*dtmin/(m_e*c**2)).to(u.Unit(''))+block.E_i).value
        
        else:
            min_index = t.index(t_min)
            u_rad = u_rad[min_index]
            E_m = (((-4/3*sigma_T*c*u_rad*block.E_i**2)*dtmin/(m_e*c**2)).to(u.Unit(''))+block.E_i).value

        block.final_energy(E_m)

    def EnergyComputer(self, blocklist, dtmin, targetlist = None, cmb = False):
        """
        Computes the energy for each block.
        
        Parameters:
        --
        blocklist :class:`~list[objects.block]`: list of blocks.
        dtmin :class:`~astropy.unit.Quantity`: time interval between the current file simulation's time and the following one.  
        targetlist :class:`~list[agnpy.target]`: list of targets that interact with the block via external compton
        cmb :class:`~boolean`: check whether to consider the CMB in the External Compton computation or not
        """

        if self.n_e == 'PowerLaw':
            for block in blocklist:
                self.E_initial(block)
                self.E_final(block, dtmin, targetlist, cmb)
                self.ElectronDistribution_PL(self.K_eNormalizer(block), block)

        elif self.n_e == 'BrokenPowerLaw':
            for block in blocklist:
                self.E_initial(block)
                self.E_final(block, dtmin, targetlist, cmb)
                self.ElectronDistribution_BPL(self.K_eNormalizer(block), block)            
        
    def K_eNormalizer(self, block):
        """
        Normalizes the K_e emission for each blob.

            K_e = rho/(integral_{E_min}^{E_final} N(E)dE * mu_0 m_p)
        """

        integral = (block.E_f**(-self.p+1)-self.gamma_min**(-self.p+1))/(-self.p+1)
        if integral == 0:
            K_e = 0*u.Unit('cm-3')
        else:
            #K_e = (block.dens/(integral*self.mu_0*(m_p.to(u.g)))).to('cm-3')
            K_e = (block.dens/(integral*self.mu_0*(m_e.to(u.g)))).to('cm-3')
        return K_e


    def ElectronDistribution_PL(self, K_e, block):
        """
        Computes the electron distribution function n_e in case of a Power Law d.f..
        Returns n_e and n_e_value for n_e((E_f+gamma_min)/2)
        """

        nelectrons = PowerLaw(
            k=K_e,
            p = self.p,
            gamma_min = self.gamma_min,
            gamma_max = block.E_f,
            mass = m_e
        )
        n_e_values = nelectrons((block.E_f+self.gamma_min)/2)
        block.electronic_distribution(K_e, nelectrons, n_e_values)

    def ElectronDistribution_BPL(self, K_e, block):
        """
        Computes the electron distribution function n_e in case of a Broken Power Law d.f..
        Returns n_e and n_e_value for n_e((E_f+gamma_min)/2)        
        """
        
        nelectrons = BrokenPowerLaw(
            k = K_e,
            p1 = self.p1,
            p2 = self.p2,
            gamma_b = self.gamma_b,
            gamma_min = self.gamma_min,
            gamma_max = block.E_f,
            mass = m_e
        )
        n_e_values = nelectrons((block.E_f+self.gamma_min)/2)
        block.electronic_distribution(K_e, nelectrons, n_e_values)