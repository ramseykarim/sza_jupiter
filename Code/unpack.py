"""
    Unpack Karto's COPSS data and the JPL Horizons distances
"""
import numpy as np
import channel as ch
# import stats as st
import subprocess
import sys

PATH = "/home/ramsey/Documents/Research/Jupiter/SZA/"
COPSS_PATH = PATH + "TXTs/jupiter_fluxes_flagged_USETHIS.txt"
HORIZ_MARS_PATH = PATH + "horizons_kartos/mars.txt"
HORIZ_JUPT_PATH = PATH + "horizons_kartos/jupt.txt"


def horizons_genfromtxt(file_name):
    """
    Columns:
        0) Julian Date
        1) RA
        2) Dec
        3) Az
        4) El
        5) LST
        6) Light Time (min)
    :param file_name: Name of ephemera file
    :return: Array of info, described above
    """
    return np.genfromtxt(file_name, delimiter=',', skip_header=63, skip_footer=59,
                         usecols=[1, 4, 5, 6, 7, 8, 9])


# noinspection SpellCheckingInspection
def unpack_horizons_new(jup, mars):
    jup_data = horizons_genfromtxt(jup)
    mars_data = horizons_genfromtxt(mars)
    return jup_data, mars_data


def unpack_copss(csv_file, frequency_list_ghz, c_type=ch.Channel):
    raw_all = np.genfromtxt(csv_file)
    date_array = raw_all[:, 0]
    channel_list = []
    for i in range(15):
        channel_list.append(c_type(raw_all[:, i * 2 + 1:i * 2 + 3], frequency_list_ghz[i]))
    return date_array, channel_list


class Unpack:
    def __init__(self):
        self.mjd = 2400000.5
        self.frequency_list_ghz = [34.688000, 34.188000, 33.688000,
                                   33.188000, 32.688000, 32.188000,
                                   31.688000, 31.188000, 30.688000,
                                   30.188000, 29.688000, 29.188000,
                                   28.688000, 28.188000, 27.688000]
        self.jupiter, self.mars = None, None
        self.dates_copss_array, self.channel_obj_list = None, None
        self.correlated_dates = None
        self.error_investigation = None
        self.channel_type = ch.Channel

    def prepare(self):
        self.jupiter, self.mars = unpack_horizons_new(HORIZ_JUPT_PATH, HORIZ_MARS_PATH)
        self.dates_copss_array, self.channel_obj_list = \
            unpack_copss(PATH + "TXTs/jupiter_fluxes_flagged_USETHIS.txt",
                         self.frequency_list_ghz, c_type=self.channel_type)
        self.correlated_dates = self.correlate_dates()
        return self

    def get_temperatures(self):
        return np.array([channel.tb for channel in self.channel_obj_list]), \
               np.array([channel.tb_error_slope for channel in self.channel_obj_list]), \
               np.array([channel.tb_error_offset for channel in self.channel_obj_list])

    def get_frequencies(self):
        return np.array(self.frequency_list_ghz)

    def get_wavelengths(self):
        return np.array([c.wavelength_cm for c in self.channel_obj_list])

    def get_distances(self):
        indices = self.correlated_dates
        c = 0.1202  # in AU/min
        distances = self.jupiter[indices, 6] * c
        return distances

    def get_dates(self):
        return self.dates_copss_array

    def get_j_el(self):
        return self.jupiter[self.correlated_dates, 4]

    def get_m_el(self):
        return self.mars[self.correlated_dates, 4]

    def correlate_dates(self):
        horizons_dates = self.jupiter[:, 0] - self.mjd
        indices = [np.argmin(np.abs(horizons_dates - x)) for x in self.dates_copss_array]
        return indices

    def adjust(self):
        self.distance_adj()
        self.alt_adj()
        self.prelim_adjust()
        self.average_adj()
        self.error_adj()
        self.synchrotron_adj()
        self.cmb_unit_adj()
        return self

    def prelim_adjust(self):
        return self

    def distance_adj(self):
        corresponding_distances = self.get_distances()
        mod_factor = (corresponding_distances / 4.04) ** 2.
        for channel in self.channel_obj_list:
            channel.distance_adj(mod_factor)
        return self

    def alt_adj(self):
        for channel in self.channel_obj_list:
            channel.altitude_adj(self.get_j_el())
        return self

    def average_adj(self):
        for channel in self.channel_obj_list:
            channel.average()
        return self

    def error_adj(self):
        flux_grid = np.array([channel.flux for channel in self.channel_obj_list])
        flux_error_grid = np.array([channel.error for channel in self.channel_obj_list])
        flux_avg_vector = np.array([channel.flux_avg for channel in self.channel_obj_list])
        day_ax = 0
        frq_ax = 1
        fractional_error_c = np.nanmean(flux_error_grid / flux_grid, axis=frq_ax) * flux_avg_vector
        flux_avg_vector = flux_avg_vector.reshape((15, 1))
        flux_grid /= flux_avg_vector
        flux_grid -= np.nanmean(flux_grid, axis=day_ax)
        flux_grid *= flux_avg_vector
        spread_daily_deviation_c = np.nanstd(flux_grid, axis=frq_ax)
        mean_daily_deviation_c = np.nanmean(flux_grid, axis=frq_ax)
        relative_error = np.sqrt(fractional_error_c**2 + spread_daily_deviation_c**2 + mean_daily_deviation_c**2)
        for channel, e in zip(self.channel_obj_list, relative_error):
            channel.error_avg = e
        return self

    def synchrotron_adj(self):
        f0 = 28.5
        j0 = 1.5
        synchrotron = j0 * ((np.array(self.frequency_list_ghz) / f0) ** (-0.4))
        for i, channel in enumerate(self.channel_obj_list):
            channel.synchrotron_adj(synchrotron[i])
        return self

    def cmb_unit_adj(self):
        for channel in self.channel_obj_list:
            channel.cmb_unit_adj()
        return self

    def error_investigation_init(self):
        self.distance_adj()
        # self.error_investigation = st.Stats(self.channel_obj_list)
        return self

    def write_raw(self, file_name):
        fl = open(file_name, 'w')
        fl.write("# Fluxes given in Jy, freq in GHz, date in MJD\n# Freq")
        for channel in self.channel_obj_list:
            msg = " %6.3f" % channel.frequency_ghz
            fl.write(msg)
        fl.write("\n# Date")
        for channel in self.channel_obj_list:
            fl.write(" Meas Err")
        for i, date in enumerate(self.dates_copss_array):
            msg = "\n%12.6f" % date
            fl.write(msg)
            for channel in self.channel_obj_list:
                msg = " %10.6f %8.6f" % (channel.flux[i], channel.error[i])
                fl.write(msg)
        fl.close()

    def write_points(self):
        fl = open('ramsey_data_07_29_17.txt', 'w')
        fl.write("# Frequency (GHz), Wavelength (cm), T_b (K), Absolute Uncertainty (K), Relative Uncertainty (K)\n")
        for f, w, t, er, ea in [channel.info_tuple() for channel in self.channel_obj_list]:
            fl.write(str(f) + ', ' + str(w) + ', ' + str(t) + ', ' + str(ea) + ', ' + str(er) + '\n')

    def print_points(self):
        for f, w, t, e, en in [channel.info_tuple() for channel in self.channel_obj_list]:
            print "F: ", f, \
                " WL: ", w, \
                "\n\tTb    : ", t, \
                "\n\tTh Er : ", e, ',', 100. * e / t, "%", \
                "\n\tEns Er: ", en, ',', 100. * en / t, "%"


def generate_names(path):
    bash_cmd = "ls " + path
    process = subprocess.Popen(["bash", "-c", bash_cmd],
                               stdout=subprocess.PIPE)
    output, error = process.communicate()
    file_name_list = output.split("\n")
    file_name_list.pop()
    return file_name_list
