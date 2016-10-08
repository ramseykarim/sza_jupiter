import matplotlib.pyplot as plt
import models as mdl
import scipy.constants as cst
import numpy as np


class Plotting:
    def __init__(self, wl=False):
        self.wl = wl
        self.labels_list = []
        self.labels_list_models = []
        plt.xlabel("Frequency (GHz)" if not self.wl else "Wavelength (cm)")
        plt.ylabel("$T_{b}$ (K)")
        plt.title("$T_{b}$ vs. " + "$\\nu$" if not self.wl else "$\lambda$")
        plt.xlim(20, 95)
        plt.ylim(120, 220)

    def add_my_data(self, unpacker):
        self.labels_list.append("SZA")
        tb = [channel.tb for channel in unpacker.channel_obj_list]
        tb_err = [channel.tb_error for channel in unpacker.channel_obj_list]
        x_axis = unpacker.frequency_list_ghz if not self.wl else [c.wavelength_cm for c in unpacker.channel_obj_list]
        plt.errorbar(x_axis, tb, yerr=tb_err, fmt='.')

    def add_model_plot(self, root_name, color, line_sty='-'):
        frequencies, points = mdl.generate_model(root_name)
        self.labels_list_models.append(root_name)
        plt.plot(self.f_or_w(frequencies), points, line_sty, color=color)

    def add_model_compare(self, root_name, color):
        frequencies, points = mdl.compare_nh3_model(root_name)
        self.labels_list_models.append((root_name + " : NH3 (new - old)"))
        plt.plot(self.f_or_w(frequencies), points, color=color)

    def add_all_other_data(self):
        for p in mdl.EXISTING_DATA:
            f, tb, tbe = p(self)
            plt.errorbar(self.f_or_w(f), tb, yerr=tbe, fmt='.')

    def f_or_w(self, frequencies):
        if self.wl:
            return (cst.c / frequencies) * 1e-7
        else:
            return frequencies

    def show_prep(self):
        plt.xscale('log', basex=2)
        tick_locations = np.arange(76) + 19
        plt.xticks(tick_locations, [str(loc) if loc % 10 == 0 else '' for loc in tick_locations])
        plt.legend(self.labels_list_models + self.labels_list, loc='lower right')

    def show(self):
        plt.show()
