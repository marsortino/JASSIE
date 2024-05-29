import numpy as np
import astropy.units as u
from astropy.constants import c, m_e, h, e, alpha, m_p, sigma_T, mu0
from agnpy.synchrotron import Synchrotron
from agnpy.compton import SynchrotronSelfCompton,ExternalCompton
from agnpy.targets import CMB
from agnpy.utils.conversion import nu_to_epsilon_prime
from agnpy.utils.math import axes_reshaper
from agnpy.compton.kernels import isotropic_kernel

# Order of classes:
# - Blumenthal
# - Bremsstrahlung
# - synchrotron
# - ExtCompton

class Bremsstrahlung:
    """

    It computes the emission produced by non thermal brehmstralung as computed by Blumenthal and Gould (1970).
    Equation (3.59):

         dN_tot/(dt*dk) = \alpha*r_0^2*c*K_e*k^{-1}* \sum_s(n_s*I_s(k, E_L;p))

    where:
         I_s(weak shielding}) \approx I'* \bar{\phi}_{weak}

    with: 
        I' = \frac{4}{3} \frac{E_L^{-(p-1)}}{p-1} - \frac{4}{3} \frac{k E_L^{-p}}{p} + \frac{k^2E_L^{-(p+1)}}{p+1}
    
    and 
                        ...4(Z^2+Z_{el})\{\ln{[2E_0(E_0-k)/k]}-1/2\}   for: k < E_0 
        \phi_weak = ...
                        ...4(Z^2+Z_{el})[\ln(4k)-1/2]                  for: k>E_0


    The emission of a relativistic flux for an isotopically
    radiating plasma blob is then computed using Dermer & Menon (2009) (5.46):

            f_nu = delta_D^4 * V'_b/(4\pi * d_l^2) * nu'  * j'_nu

    primed quantities are in the comoving spaceframe of the blob.
    nu' is in units of m_e*c^2
    """

    def __init__(self, blocklist, config):
        """
        Parameters:
        --
        blocklist :class:`~list[objects.block]`: list of blocks
        nu :class:`~astropy.quantity`: numpy array of frequency
        """
        
        nu = config['nu']
        self.gamma_min = config['gamma_Min']
        self.B = config['B']
        self.p = config['p']
        self.Z = config['Z']
        self.r_e = (e.esu**2/(m_e*c**2)).to(u.m)
        self.k = (h*nu/(m_e*c**2)).to(u.Unit(""))
        Total_J_nu = self.emission(blocklist, nu)
        sed = self.computeflux(blocklist, nu, Total_J_nu)
        self.sedlist = sed[0]
        self.total_sed = sed[1]


    def emission(self, blocklist, nu):
        """
        Computes the bremsstrahlung emission. Returns: 
                alpha * r_e^2*c*K_e*k^-1 * integral
        integral is computed in the function 'integral'

        Parameters:
        --
        blocklist :class:`~list[objects.block]`: list of blocks
        nu :class:`~astropy.quantity`: numpy array of frequency
        """


        TotalIntegral = self.integral()
        TotalIntegral = np.where(nu>1e13*u.Hz, 0, TotalIntegral)

        BlocksEmission = []

        for block in blocklist:
            prefactor = alpha*self.r_e**2
            midfactor = block.K_e*c
            Integral = self.k**(-1)*block.n_0*TotalIntegral
            rest_energy = m_e*c**2
            total = prefactor*midfactor*Integral*rest_energy
            total = total.to(u.Unit('erg cm-3 s-1'))
            BlocksEmission.append(total)
        return BlocksEmission

    def integral(self):
        """
        Computes the I' value of Blumenthal & Gould (1970) equation "I = I'*\phi" as it follows:
            I' = \frac{4}{3} \frac{E_L^{-(p-1)}}{p-1} - \frac{4}{3} \frac{k E_L^{-p}}{p} + \frac{k^2E_L^{-(p+1)}}{p+1}
        """
        p = self.p
        k = self.k

        E_l = np.maximum.outer(k, self.gamma_min)
        factor = (4/3*E_l**(-(p-1))/(p-1) - 4/3*k*E_l**(-p)/p + k**2*E_l**(-p-1)/(p+1))*self.phi_weak(k)
        return factor
        
    def phi_weak(self, k):
        """ 
        From Blumenthanl and Gould (1970), equations [3.64] for phi_weak:
                        ...4(Z^2+Z_{el})\{\ln{[2E_0(E_0-k)/k]}-1/2\}   for: k < E_0 
        \phi_weak = ...
                        ...4(Z^2+Z_{el})[\ln(4k)-1/2]                  for: k >= E_0
        Note that: for the purpose of the computation Z_el = Z
        
        """
        E_0 = self.gamma_min
        Z = self.Z
        return np.where(k>=E_0, 4*(Z**2+Z)*(np.log(4*k)-1/2),  4*(Z**2+Z)*(np.log(2*E_0*(E_0-k)/k)-0.5))


    def computeflux(self, blockslist, nu, Total_J_nu):
        """
        Computes the flux for each blob in blobslist

        Parameters:
        --
        blocklist :class:`~list[objects.block]`: list of blocks
        nu :class:`~astropy.quantity`: numpy array of frequency

        Total_J_nu :class:`~emission.Bremmstrahlung.emission`: non-thermal bremmstrahlung emission computed for each block.
        """
        sedlist = []
        totalsed = 0
        for block, J_nu in zip(blockslist, Total_J_nu):
            sed = self.flux(block, nu, J_nu)
            sedlist.append(sed)
            totalsed +=sed
        return sedlist, totalsed

    def flux(self, block, nu, J_nu):
        """
        Computes the flux using Dermer (5.46) for an isotropic radiation from a relativistic plasma blob:
                    f_nu = delta_D^4 * V'_b/(4\pi * d_l^2) * nu'  * j'_nu


        """
        V_b = 4/3*np.pi*block.radius**3
        epsilon = h*nu/(m_e*c**2)
        return (block.delta_D**4*V_b/(4*np.pi*block.obs_distance**2)*epsilon*J_nu).to(u.Unit('erg cm-2 s-1'))
    
class synchrotron:
    """
    Computes the SED produced by each Blob using AGNpy routines.
    """
    def __init__(self, blocklist, config, whole_blocklist):
        """
        Parameters:
        --
        blocklist :class:`~list[objects.block]`: list of blocks
        config :class:`~dictionary`: Dictionary containing info on B and distribution function parameters.
        """
        self.B = config['B']
        nu = config['nu']

        if config['n_e'] == 'PowerLaw':
            self.p = config['p']
        elif config['n_e'] == 'BrokenPowerLaw': # nb. broken power law is still not available.
            self.p1 = config['p1']
            self.p2 = config['p2']
            self.gamma_B = config['gamma_B']

        self.ssa = config['id_ssa']
        self.ssc = config['id_ssc']
        self.whole_blocklist = whole_blocklist
        id_ssc = config['id_ssc']
        id_ssa = config['id_ssa']

        if id_ssc == False:        
            if id_ssa:
                self.flux_ssa(blocklist, nu)
            else:
                self.flux(blocklist, nu)

        elif id_ssc == True:
            if id_ssa:
                self.flux_ssc_ssa(blocklist, nu)
            else:
                self.flux_ssc(blocklist, nu)


    def flux(self, blocklist, nu):
        """
        Computes the flux in the case in which the self compton is not present.
        """
        synchro = 0

        for block in blocklist:
            sed = Synchrotron(block.blob).sed_flux(nu)
            synchro += sed*np.exp(-block.k*block.obs_raypath)
        self.total_sed = synchro

    def flux_ssa(self, blocklist, nu):

        synchro_ssa = 0
        for block in blocklist:
            sed = Synchrotron(block.blob, ssa = True).sed_flux(nu)
            attenuation = 0
            for index in block.distance_order:
                absorbing_blob = blocklist[index].blob
                attenuation += self.synchro_self_absorption(nu, absorbing_blob)
            synchro_ssa += sed*np.exp(-block.k*block.obs_raypath)*attenuation

        self.total_sed = synchro_ssa

    def flux_ssc(self, blocklist, nu):
        """
        Computes the flux in the case in which the self compton is present.
        """
        synchro = 0
        synchro_ssc = 0
        for block in blocklist:
            sed = Synchrotron(block.blob).sed_flux(nu)
            synchro += sed*np.exp(-block.k*block.obs_raypath)
            ssc = SynchrotronSelfCompton(block.blob).sed_flux(nu)
            ssc_by_blobs_on_lov = 0
            for index in block.distance_order:
                scattering_blob = self.whole_blocklist[index].blob
                ssc_sed = self.synchroselfcompton(nu, block.blob, scattering_blob, *scattering_blob.n_e.parameters)
                ssc_by_blobs_on_lov += self.ray_path(ssc_sed, self.whole_blocklist[index])
            synchro_ssc += self.ray_path(ssc, block) + ssc_by_blobs_on_lov
            
        self.total_sed = synchro+synchro_ssc
    
    def flux_ssc_ssa(self, blocklist, nu):
        """
        Computes the flux in the case in which the self compton is present.
        """
        synchro_ssa = 0
        synchro_ssc = 0
        print('Starting computing Synchrotron and SSC, it may take some time.')
        for block in blocklist:
            sed = Synchrotron(block.blob, ssa = True).sed_flux(nu)

            ssc = SynchrotronSelfCompton(block.blob, ssa = True).sed_flux(nu)
            synchro_ssc += self.ray_path(ssc, block)
            ssc_by_blobs_on_lov = 0
            attenuation = 0
            for index in block.distance_order:
                scattering_blob = self.whole_blocklist[index].blob
                ssc_sed = self.synchroselfcompton(nu, block.blob, scattering_blob, *scattering_blob.n_e.parameters)
                attenuation += self.synchro_self_absorption(nu, scattering_blob, *scattering_blob.n_e.parameters)
                ssc_by_blobs_on_lov += self.ray_path(ssc_sed, self.whole_blocklist[index])
            synchro_ssa += sed*np.exp(-block.k*block.obs_raypath)*attenuation
    
        self.total_sed = synchro_ssa+synchro_ssc

    def ray_path(self, sed, block):
        """
        Returns the specific intensity after that it went through the blobs on the line of view.
        """
        return sed*np.exp((-block.k*block.obs_raypath).to(u.Unit('')))
    
    def synchroselfcompton(self, nu, emitting_blob, scattering_blob, *args, integrator = np.trapz,):
        """
        Computes the synchrotron self Compton, using agnpy method and routines. Modified in order to take
        into account the synchrotron emitted by a different blob.
        """
        nu_to_integrate = np.logspace(5, 30, 200) * u.Hz  # used for SSC
        epsilon = nu_to_epsilon_prime(nu_to_integrate, emitting_blob.z, emitting_blob.delta_D)
        epsilon_s = nu_to_epsilon_prime(nu, scattering_blob.z, scattering_blob.delta_D)

        emitting_blob_synchro = Synchrotron(emitting_blob).evaluate_sed_flux(nu_to_integrate,
                                                                             emitting_blob.z,
                                                                             emitting_blob.d_L,
                                                                             emitting_blob.delta_D,
                                                                             emitting_blob.B,
                                                                             emitting_blob.R_b,
                                                                             emitting_blob.n_e,
                                                                             *args,                                                                             
                                                                             ssa = self.ssa,
                                                                             integrator = integrator,
                                                                             gamma = emitting_blob.gamma_e,
        )

        ## check u_synch eq.
        u_synch = (3 * np.power(scattering_blob.d_L, 2) * emitting_blob_synchro) / (
            c * np.power(scattering_blob.R_b, 2) * np.power(scattering_blob.delta_D, 4) * epsilon
        )
        _gamma, _epsilon, _epsilon_s = axes_reshaper(scattering_blob.gamma_e, epsilon, epsilon_s)
        V_b = 4 / 3 * np.pi * np.power(scattering_blob.R_b, 3)
        N_e = V_b * scattering_blob.n_e.evaluate(_gamma, *args)
        # reshape u as epsilon
        _u_synch = np.reshape(u_synch, (1, u_synch.size, 1))
        # integrate
        kernel = isotropic_kernel(_gamma, _epsilon, _epsilon_s)
        integrand = (
            _u_synch / np.power(_epsilon, 2) * N_e / np.power(_gamma, 2) * kernel
        )
        integral_gamma = integrator(integrand, scattering_blob.gamma_e, axis=0)
        integral_epsilon = integrator(integral_gamma, epsilon, axis=0)
        emissivity = 3 / 4 * c * sigma_T * np.power(epsilon_s, 2) * integral_epsilon
        prefactor = np.power(scattering_blob.delta_D, 4) / (4 * np.pi * np.power(scattering_blob.d_L, 2))
        return (prefactor * emissivity).to("erg cm-2 s-1")

    def synchro_self_absorption(self, nu, absorbing_blob, *args):
        """
        Computes the self absorption for the blobs standing between the emitting blob and the observer using agnpy method.
        """
        tau = Synchrotron.evaluate_tau_ssa(nu, 
                                    absorbing_blob.z,
                                    absorbing_blob.d_L, 
                                    absorbing_blob.delta_D, 
                                    absorbing_blob.B, 
                                    absorbing_blob.R_b,
                                    absorbing_blob.n_e,
                                    *args,
                                    gamma = absorbing_blob.gamma_e,
                                    )
        attenuation = tau_to_attenuation(tau)
        return attenuation
    
class ExtCompton:
    """
    Computes the external compton scattering for each blob using AGNpy.
    """
    def __init__(self, blocklist, nu, targetlist = [], id_cmb = False):
        """
        Parameters:
        --
        blocklist :class:`~list[objects.block]`: list of blocks
        nu :class:`~astropy.quantity`: Array of frequencies over which the SED needs to be computed.
        targetlist :class:`~list[agnpy.targets]`:list of target over which the external compton is computed.
        id_cmb :class:`~boolean`: True/False whether to consider the CMB or not.
        """

        if id_cmb == False:
            self.total_sed = self.totalflux_no_cmb(blocklist, nu, targetlist)
        elif id_cmb == True:
            self.total_sed = self.totalflux(blocklist, nu, targetlist)

    def totalflux(self, blocklist, nu, targetlist):
        """
        Computes the EC emission - CMB included - over each blob and returns the total sed. 
        """
        totalsed = 0
        i_ec = 0
        i_cmb = 0

        for block in blocklist:
            sed_list = []

            if block.EC:
                i_ec += 1
                for target in targetlist:
                    sed_list.append(ExternalCompton(block.blob, target, block.distance))
                    
            if block.CMB:
                i_cmb += 1
                ec = ExternalCompton(block.blob, CMB(block.redshift), block.distance)
                totalsed +=  self.lineofview(ec.sed_flux(nu), block)

            for sed in sed_list:
                sed1 = sed.sed_flux(nu)
                totalsed += self.lineofview(sed1, block)

        print('total CMB SED computed:', i_cmb)
        print('total EC SED computed:', i_ec)
        return totalsed
    
    def totalflux_no_cmb(self, blocklist,nu, targetlist):
        """
        Computes the EC emission - NO CMB - over each blob and returns the total sed. 
        """
        totalsed = 0
        i = 0
        for block in blocklist:
            
            sed_list = []
        
            if block.EC:
                i += 1
                for target in targetlist:
                    sed_list.append(ExternalCompton(block.blob, target, block.distance))

            for sed in sed_list:
                sed1 = sed.sed_flux(nu)
                
                totalsed +=  self.lineofview(sed1, block)

        print('total EC SED computed:', i)
        return totalsed
    
    def lineofview(self, sed, block):
        return sed*np.exp((-block.k*block.obs_raypath).to(u.Unit('')))

def tau_to_attenuation(tau):
    """FROM AGNPY ***
    Converts the synchrotron self-absorption optical depth to an attenuation
    Eq. 7.122 in [DermerMenon2009]_."""
    u = 1 / 2 + np.exp(-tau) / tau - (1 - np.exp(-tau)) / np.power(tau, 2)
    return np.where(tau < 1e-3, 1, 3 * u / tau)
