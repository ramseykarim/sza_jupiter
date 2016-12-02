"""
    Class prototype for a channel
"""
import scipy.constants as cst
import numpy as np
import copy

jy_to_w = 1.e-26
c = cst.c  # m/s
h = cst.h  # Planck, SI
pi = cst.pi
k = cst.Boltzmann  # J/K
T_cmb = 2.725  # Kelvins


class Channel:

    def __init__(self, rectangle, frequency_ghz):
        self.flux = rectangle[:, 0]
        self.error = rectangle[:, 1]  # Original error from txt, assumed in janskys
        self.frequency_ghz = frequency_ghz
        self.wavelength_cm = (cst.c / frequency_ghz) * 1e-7
        self.flux_avg = 0.
        self.error_avg = 0.  # Janskys, very small, good for slope
        self.ensemble_error = 0.  # Janskys, ~5%, good for absolute calibration
        self.tb = 0.
        self.tb_error_slope = 0.
        self.tb_error_offset = 0.

    def distance_adj(self, mod_array):
        self.flux *= mod_array

    def average(self):
        """
        This is the newest average, corresponding to:
        Weight:
        w = 1/error^2
        Flux:
        <F> = (sum F * w) / (sum w)
        Error:
        sigma_T^2 = sigma^2 + n^2
        where sigma is
        sigma^2 = (sum (F - <F>)^2 * w) / (sum w)
        and n is
        n^2 = 1/(sum w)
        :return: The averages
        """
        assert isinstance(self.error, np.ndarray)
        assert isinstance(self.flux, np.ndarray)
        # weights = 1. / (self.error * self.flux) ** 2.  # This is the original
        weights = 1. / self.error ** 2.  # This is new
        self.flux_avg = np.nansum(weights * self.flux) / np.nansum(weights)
        variance_sq = 1. / np.nansum(weights)
        accuracy_error_sq = np.nansum(weights * (self.flux - self.flux_avg)**2.) * variance_sq
        self.ensemble_error = np.sqrt(accuracy_error_sq + variance_sq)
        # total_error_sq_inv = (1. / accuracy_error_sq) + (1. / variance_sq)
        # self.error_avg = 1. / np.sqrt(total_error_sq_inv)
        self.error_avg = np.nanmean(self.error/self.flux) * self.flux_avg

    def synchrotron_adj(self, synchrotron):
        self.flux_avg -= synchrotron

    def cmb_unit_adj(self):
        f_i = self.flux_avg * jy_to_w
        sigma_i = self.error_avg * jy_to_w
        sigma_hat = self.ensemble_error * jy_to_w
        v = self.frequency_ghz * 1.e9
        r_eq = 71492000.
        r_p = 66854000.
        sub_earth_lat = 0.15 * cst.degree  # radians
        r_pp = np.sqrt((r_eq * np.sin(sub_earth_lat)) ** 2. + (r_p * np.cos(sub_earth_lat)) ** 2.)
        first_term = 1. / (np.exp((h * v) / (k * T_cmb)) - 1.)
        second_term_cst = (c ** 2. * (cst.au * 4.04)**2.) / (2. * h * v**3. * pi * r_eq * r_pp)

        second_term_flx = f_i * second_term_cst
        second_term_hi_err = (f_i + sigma_i) * second_term_cst
        second_term_lo_err = (f_i - sigma_i) * second_term_cst
        second_term_hi_ens_err = (f_i + sigma_hat) * second_term_cst
        second_term_lo_ens_err = (f_i - sigma_hat) * second_term_cst

        inside_log_flx = 1. + 1. / (first_term + second_term_flx)
        inside_log_hi_err = 1. + 1. / (first_term + second_term_hi_err)
        inside_log_lo_err = 1. + 1. / (first_term + second_term_lo_err)
        inside_log_hi_ens_err = 1. + 1. / (first_term + second_term_hi_ens_err)
        inside_log_lo_ens_err = 1. + 1. / (first_term + second_term_lo_ens_err)

        self.tb = ((h * v) / k) / np.log(inside_log_flx)
        tb_hi_error = ((h * v) / k) / np.log(inside_log_hi_err)
        tb_lo_error = ((h * v) / k) / np.log(inside_log_lo_err)
        hi_error, lo_error = np.abs(tb_hi_error - self.tb), np.abs(tb_lo_error - self.tb)
        tb_hi_ens_error = ((h * v) / k) / np.log(inside_log_hi_ens_err)
        tb_lo_ens_error = ((h * v) / k) / np.log(inside_log_lo_ens_err)
        hi_ens_error, lo_ens_error = np.abs(tb_hi_ens_error - self.tb), np.abs(tb_lo_ens_error - self.tb)
        self.tb_error_slope = np.max([hi_error, lo_error])  # This was multiplied by 3 originally? #TODO
        self.tb_error_offset = np.max([hi_ens_error, lo_ens_error])

    def altitude_adj(self, original_alt):
        finite_i = np.isfinite(self.flux)
        flux = self.flux[finite_i]
        alt = original_alt[finite_i]
        degree = 1
        alt_fit = np.polyfit(alt, flux, deg=degree)
        alt_solution = np.zeros(len(alt))
        for i, coefficient in enumerate(alt_fit[:degree]):
            alt_solution += coefficient * (alt ** (degree - i))
        self.flux[finite_i] -= alt_solution

    def info_tuple(self):
        return (self.frequency_ghz,
                self.wavelength_cm,
                self.tb,
                self.tb_error_slope,
                self.tb_error_offset)

    def copy(self):
        return copy.deepcopy(self)

