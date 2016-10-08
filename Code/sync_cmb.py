import numpy as np
import scipy.constants as cst


jy_to_w = 1e-26
c = cst.c  # m/s
h = cst.h  # Planck, SI
pi = cst.pi
k = cst.Boltzmann  # J/K
T_cmb = 2.725  # Kelvins


def synchrotron_cmb(point, frequency):
    new_p = point
    f0 = 28.5
    j0 = 1.5
    synchrotron = j0 * ((frequency / f0) ** (-0.4))
    new_p -= synchrotron
    f_i = new_p * jy_to_w
    v = frequency * 1e9
    r_eq = 71492000.
    r_p = 66854000.
    sub_earth_lat = 0.15 * cst.degree  # radians
    r_pp = np.sqrt((r_eq * np.sin(sub_earth_lat)) ** 2. + (r_p * np.cos(sub_earth_lat)) ** 2.)
    first_term = 1. / (np.exp((h * v) / (k * T_cmb)) - 1)
    second_term_cst = (c ** 2. * (cst.au * 4.04) ** 2.) / (2 * h * v ** 3. * pi * r_eq * r_pp)
    second_term_flx = f_i * second_term_cst
    inside_log_flx = 1. + 1. / (first_term + second_term_flx)
    tb = ((h * v) / k) / np.log(inside_log_flx)
    return tb


synchrotron_cmb_units_np = np.vectorize(synchrotron_cmb)
