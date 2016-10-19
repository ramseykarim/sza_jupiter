"""
    Unpack Karto's COPSS data and the JPL Horizons distances
"""
import numpy as np
import csv
import channel as ch
import matplotlib.pyplot as plt
import stats as st

PATH = "/home/ramsey/Documents/Research/Jupiter/SZA/"


def unpack_horizons(csv_file):
    data_file = open(csv_file, 'rb')
    reader = csv.reader(data_file, delimiter=',')
    distance_list = []
    date_list = []
    for row in reader:
        distance_list.append(float(row[len(row) - 2]))
        date_list.append(float(row[1]))
    c = 0.1202  # in AU/min
    distance_array = np.array(distance_list)
    distance_array *= c
    date_array = np.array(date_list)
    return date_array, distance_array


def unpack_copss(csv_file, frequency_list_ghz):
    raw_all = np.genfromtxt(csv_file)
    date_array = raw_all[:, 0]
    channel_list = []
    for i in range(15):
        channel_list.append(ch.Channel(raw_all[:, i*2 + 1:i*2 + 3], frequency_list_ghz[i]))
    return date_array, channel_list


class Unpack:

    def __init__(self):
        self.mjd = 2400000.5
        self.frequency_list_ghz = [34.688000, 34.188000, 33.688000,
                                   33.188000, 32.688000, 32.188000,
                                   31.688000, 31.188000, 30.688000,
                                   30.188000, 29.688000, 29.188000,
                                   28.688000, 28.188000, 27.688000]
        # Unpack data files into object
        self.dates_horizons_array, self.distances_horizons_array =\
            unpack_horizons(PATH + "horizons_kartos/distances.csv")
        self.dates_horizons_array -= self.mjd
        self.dates_copss_array, self.channel_obj_list =\
            unpack_copss(PATH + "TXTs/jupiter_fluxes_flagged_USETHIS.txt",
                         self.frequency_list_ghz)
        self.error_investigation = 0

    def adjust(self):
        self.distance_adj()
        self.average_adj()
        self.synchrotron_adj()
        self.cmb_unit_adj()
        return self

    def distance_adj(self):
        indices = [np.argmin(np.abs(self.dates_horizons_array - x)) for x in self.dates_copss_array]
        corresponding_distances = np.array([self.distances_horizons_array[i] for i in indices])
        mod_factor = (corresponding_distances / 4.04)**2.
        for channel in self.channel_obj_list:
            channel.distance_adj(mod_factor)

    def average_adj(self):
        for channel in self.channel_obj_list:
            channel.average()

    def synchrotron_adj(self):
        f0 = 28.5
        j0 = 1.5
        synchrotron = j0 * ((np.array(self.frequency_list_ghz) / f0)**(-0.4))
        for i, channel in enumerate(self.channel_obj_list):
            channel.synchrotron_adj(synchrotron[i])

    def cmb_unit_adj(self):
        for channel in self.channel_obj_list:
            channel.cmb_unit_adj()

    def plot_simple(self, wl=False):
        tb = [channel.tb for channel in self.channel_obj_list]
        tb_err = [channel.tb_error for channel in self.channel_obj_list]
        x_axis = self.frequency_list_ghz if not wl else [c.wavelength_cm for c in self.channel_obj_list]
        plt.errorbar(x_axis, tb, yerr=tb_err, fmt='.')
        plt.xlabel("Frequency (GHz)" if not wl else "Wavelength (cm)")
        plt.ylabel("$T_{b}$ (K)")
        plt.title("$T_{b}$, synchrotron & CMB corrected")
        plt.show()

    def error_investigation_init(self):
        self.distance_adj()
        self.error_investigation = st.Stats(self.channel_obj_list)
        return self

    def write_points(self):
        fl = open('ramsey_data_10_17_16.txt', 'w')
        fl.write("# Frequency (GHz), Wavelength (cm), T_b (K), Thermal Error (K), Ensemble Error (K)\n")
        for f, w, t, e, en in [channel.info_tuple() for channel in self.channel_obj_list]:
            fl.write(str(f) + ', ' + str(w) + ', ' + str(t) + ', ' + str(e) + ', ' + str(en) + '\n')

    def print_points(self):
        for f, w, t, e, en in [channel.info_tuple() for channel in self.channel_obj_list]:
            print "F: ", f, " WL: ", w, "\n\tTb    : ", t, "\n\tTh Er : ", e, "\n\tEns Er: ", en
