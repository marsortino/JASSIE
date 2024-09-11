import numpy as np
import astropy.units as u
from astropy.constants import c, sigma_T, G
from agnpy.utils.math import (
    axes_reshaper,
    gamma_e_to_integrate,
    phi_to_integrate,
)
from agnpy.utils.conversion import nu_to_epsilon_prime, to_R_g_units
from agnpy.targets import(
    SSDisk,
)
from agnpy.compton import compton_kernel



class ExternalCompton:
    """
    
    """
    def __init__(self, block, target, integrator = np.trapz):
        self.block = block
        self.blob = block.blob
        self.integrator = integrator
        self.target_settings(target)

    def target_settings(self, target):
        if isinstance(target, SSDisk):
            self.disk = target

    def sed_flux_from_disk(self, nu, mu_s):
        flux = self.flux_from_disk(
            nu,
            self.blob.z,
            self.blob.d_L,
            self.blob.delta_D,
            mu_s,
            self.blob.R_b,
            self.block.x*u.cm,
            self.block.y*u.cm,
            self.block.z*u.cm,
            self.disk.M_BH,
            self.disk.L_disk,
            self.disk.eta,
            self.disk.R_in,
            self.disk.R_out,
            self.blob.n_e,
            *self.blob.n_e.parameters,
            integrator=self.integrator,
            gamma=self.blob.gamma_e_external_frame
        )
        return flux

    def flux_from_disk(self,
        nu,
        z,
        d_L, 
        delta_D, 
        mu_s, 
        R_b, 
        x_blob, 
        y_blob,
        z_blob, 
        M_BH, 
        L_disk, 
        eta, 
        R_in, 
        R_out, 
        n_e, 
        *args,
        integrator = np.trapz, 
        gamma = gamma_e_to_integrate, 
        mu_size = 100, 
        ):
        """
        
        """

        epsilon_s = nu_to_epsilon_prime(nu, z)
        r_tilde = to_R_g_units(np.abs(z_blob), M_BH)
        R_out_tilde = to_R_g_units(R_out, M_BH)
        R_in_tilde = to_R_g_units(R_in, M_BH)

        m_dot = (L_disk / (eta * np.power(c, 2))).to("g /s ")

        if x_blob != 0 and y_blob != 0:
            if x_blob > 0:
                phi_1, phi_2 = splitted_phi_to_integrate()
            else:
                phi_2, phi_1 = splitted_phi_to_integrate()
            
            distance_blob_z_axis = np.sqrt(np.power(x_blob, 2) + np.power(y_blob, 2))
            distance_blob_z_axis_tilde = to_R_g_units(distance_blob_z_axis, M_BH)
            mu_1, mu_2 = self.evaluate_mu_disk(R_in_tilde, R_out_tilde, r_tilde, distance_blob_z_axis_tilde, mu_size)
            _gamma, _mu_1, _phi_1, _epsilon_s = axes_reshaper(gamma, mu_1, phi_1, epsilon_s)
            _gamma, _mu_2, _phi_2, _epsilon_s = axes_reshaper(gamma, mu_2, phi_2, epsilon_s) # Not ideal, to be changed.
            V_b = 4 / 3 * np.pi * np.power(R_b, 3)
            N_e = V_b * n_e.evaluate(_gamma / delta_D, *args) 

            epsilon_1 = SSDisk.evaluate_epsilon_mu(L_disk, M_BH, eta, _mu_1, r_tilde)

            epsilon_2 = SSDisk.evaluate_epsilon_mu(L_disk, M_BH, eta, _mu_2, r_tilde)
            
            phi_disk_1 = SSDisk.evaluate_phi_disk_mu(_mu_1, R_in_tilde, r_tilde)

            phi_disk_2 = SSDisk.evaluate_phi_disk_mu(_mu_2, R_in_tilde, r_tilde)

            kernel_1 = compton_kernel(_gamma, _epsilon_s, epsilon_1, mu_s, _mu_1, _phi_1)

            kernel_2 = compton_kernel(_gamma, _epsilon_s, epsilon_2, mu_s, _mu_2, _phi_2)

            integrand_1 = (
                phi_disk_1
                / np.power(epsilon_1, 2)
                / _mu_1
                / np.power(np.power(_mu_1, -2) - 1, 3 / 2)
                * N_e
                / np.power(_gamma, 2)
                * kernel_1
            )  
            integrand_2 = (
                phi_disk_2
                / np.power(epsilon_2, 2)
                / _mu_2
                / np.power(np.power(_mu_2, -2) - 1, 3 / 2)
                * N_e
                / np.power(_gamma, 2)
                * kernel_2
            )   
            integral_gamma_1 = integrator(integrand_1, gamma, axis=0)
            integral_mu_1 = np.trapz(integral_gamma_1, mu_1, axis=0)
            integral_phi_1 = np.trapz(integral_mu_1, phi_1, axis=0)

            integral_gamma_2 = integrator(integrand_2, gamma, axis=0)
            integral_mu_2 = np.trapz(integral_gamma_2, mu_2, axis=0)
            integral_phi_2 = np.trapz(integral_mu_2, phi_2, axis=0)
            
            integral_phi = integral_phi_1 + integral_phi_2

        else:
            phi = phi_to_integrate()
            mu = SSDisk.evaluate_mu_from_r_tilde(R_in_tilde, R_out_tilde, r_tilde, mu_size)
            _gamma, _mu, _phi, _epsilon_s = axes_reshaper(gamma, mu, phi, epsilon_s)
            V_b = 4 / 3 * np.pi * np.power(R_b, 3)
            N_e = V_b * n_e.evaluate(_gamma / delta_D, *args)
            epsilon = SSDisk.evaluate_epsilon_mu(L_disk, M_BH, eta, _mu, r_tilde)
            phi_disk = SSDisk.evaluate_phi_disk_mu(_mu, R_in_tilde, r_tilde)
            kernel = compton_kernel(_gamma, _epsilon_s, epsilon, mu_s, _mu, _phi)
            integrand = (
                phi_disk
                / np.power(epsilon, 2)
                / _mu
                / np.power(np.power(_mu, -2) - 1, 3 / 2)
                * N_e
                / np.power(_gamma, 2)
                * kernel
            )
            integral_gamma = integrator(integrand, gamma, axis=0)
            integral_mu = np.trapz(integral_gamma, mu, axis=0)
            integral_phi = np.trapz(integral_mu, phi, axis=0)
        
        prefactor_num = (
            9
            * sigma_T
            * G
            * M_BH
            * m_dot
            * np.power(epsilon_s, 2)
            * np.power(delta_D, 3)
        )

        prefactor_denom = (
            np.power(2, 9) * np.power(np.pi, 3) * np.power(d_L, 2) * np.power(np.abs(z_blob), 3)
        )
        return (prefactor_num / prefactor_denom * integral_phi).to("erg cm-2 s-1")


    def evaluate_mu_disk(self, R_in_tilde, R_out_tilde, r_tilde, distance_blob_z_tilde, size = 100):
        """
        
        """
        abs_dist = np.abs(distance_blob_z_tilde)
        mu_1_max =  1 / np.sqrt(1+np.power( ( (R_out_tilde + abs_dist) / r_tilde ), 2) )
        mu_1_min =  1 / np.sqrt(1+np.power( ( (R_in_tilde + abs_dist) / r_tilde ), 2) )
        # Checks which angle is smaller and assigns the correct variable names.
        if mu_1_max < mu_1_min:
            tmp = mu_1_max
            mu_1_max = mu_1_min
            mu_1_min = tmp

        mu_2_max =  1 / np.sqrt(1+np.power( ( (R_out_tilde - abs_dist) / r_tilde ), 2) )
        mu_2_min =  1 / np.sqrt(1+np.power( ( (R_in_tilde - abs_dist) / r_tilde ), 2) )
        # Checks which angle is smaller and assigns the correct variable names.
        if mu_2_max < mu_2_min:
            tmp = mu_2_max
            mu_2_max = mu_2_min
            mu_2_min = tmp

        return (np.linspace(mu_1_min, mu_1_max, int(size/2)), np.linspace(mu_2_min, mu_2_max, int(size/2)))


def splitted_phi_to_integrate():
    return [np.linspace(0, np.pi, 25), np.linspace(np.pi, 2*np.pi, 25)]