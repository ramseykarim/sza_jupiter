"""
    Class prototype for a channel
"""
import scipy.constants as cst
import numpy as np

jy_to_w = 1e-26
c = cst.c  # m/s
h = cst.h  # Planck, SI
pi = cst.pi
k = cst.Boltzmann  # J/K
T_cmb = 2.725  # Kelvins


class Channel:

    def __init__(self, rectangle, frequency_ghz):
        self.flux = rectangle[:, 0]
        self.error = rectangle[:, 1]  # Still a percentage
        self.frequency_ghz = frequency_ghz
        self.wavelength_cm = (cst.c / frequency_ghz) * 1e-7
        self.flux_avg = 0.
        self.error_avg = 0.  # Janskys
        self.tb = 0.
        self.tb_error = 0.

    def distance_adj(self, mod_array):
        self.flux *= mod_array

    def average(self):
        assert isinstance(self.error, np.ndarray)
        assert isinstance(self.flux, np.ndarray)

        weights = 1. / (self.error * self.flux) ** 2.
        self.flux_avg = np.nansum(weights * self.flux) / np.nansum(weights)
        variance_sq = 1. / np.nansum(weights)
        acc_error_sq = np.nansum(weights * (self.flux - self.flux_avg)**2.) * variance_sq
        total_error_sq_inv = (1. / acc_error_sq) + (1. / variance_sq)
        self.error_avg = 1. / np.sqrt(total_error_sq_inv)

    def average_old(self):
        assert isinstance(self.error, np.ndarray)
        assert isinstance(self.flux, np.ndarray)
        weights = 1. / (self.error * self.flux)
        self.flux_avg = np.nansum(1. / self.error) / np.nansum(weights)
        self.error_avg = 1. / np.sqrt(np.nansum(((self.flux - self.flux_avg) * weights)**2.))

    def synchrotron_adj(self, synchrotron):
        self.flux_avg -= synchrotron

    def cmb_unit_adj(self):
        f_i = self.flux_avg * jy_to_w
        sigma_i = self.error_avg * jy_to_w
        v = self.frequency_ghz * 1e9
        r_eq = 71492000.
        r_p = 66854000.
        sub_earth_lat = 0.15 * cst.degree  # radians
        r_pp = np.sqrt((r_eq * np.sin(sub_earth_lat)) ** 2. + (r_p * np.cos(sub_earth_lat)) ** 2.)
        first_term = 1. / (np.exp((h * v) / (k * T_cmb)) - 1)
        second_term_cst = (c ** 2. * (cst.au * 4.04)**2.) / (2 * h * v**3. * pi * r_eq * r_pp)
        second_term_flx = f_i * second_term_cst
        second_term_err = sigma_i * second_term_cst
        inside_log_flx = 1. + 1. / (first_term + second_term_flx)
        inside_log_err = 1. + 1. / second_term_err
        self.tb = ((h * v) / k) / np.log(inside_log_flx)
        self.tb_error = ((h * v) / k) / np.log(inside_log_err)


