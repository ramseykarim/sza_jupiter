import matplotlib.pyplot as plt
import models as mdl
import unpack as up
import datacontainer as dc
import scipy.constants as cst
import numpy as np


class Plotting:
    def __init__(self, wl=False):
        self.wl = wl
        self.labels_list = []
        self.labels_list_models = []

    def prep_data_plot(self):
        plt.xlabel("$\nu$ (GHz)" if not self.wl else "$\lambda$ (cm)")
        plt.ylabel("$T_{b}$ (K)")
        plt.title("$T_{b}$ vs. " + "$\\nu$" if not self.wl else "$\lambda$")
        plt.xlim(20, 95)
        plt.ylim(120, 220)

    def add_my_data(self, unpacker):
        if isinstance(unpacker, dc.DataContainer):
            x_axis, tb, tb_err, tb_ens_err = unpacker.yield_model_comp_info()
        else:
            self.labels_list.append("This Work")
            tb = [channel.tb for channel in unpacker.channel_obj_list]
            tb_err = [channel.tb_error_slope for channel in unpacker.channel_obj_list]
            tb_ens_err = [channel.tb_error_offset for channel in unpacker.channel_obj_list]
            x_axis = unpacker.frequency_list_ghz if not self.wl else [c.wavelength_cm for c in unpacker.channel_obj_list]
        plt.errorbar(x_axis, tb, yerr=tb_ens_err, fmt='o', color='black', mec='black', mfc='white', mew=1.2)
        # plt.errorbar(x_axis, tb, yerr=tb_err, fmt='.', label='_nolegend_')

    def add_model_plot(self, root_name, color, line_sty='-', label=None, linewidth=2.0):
        if label is None:
            label = root_name
        frequencies, points = mdl.generate_model(root_name)
        x_values = self.f_or_w(frequencies)
        if color.lower() == "default":
            plt.plot(x_values, points, line_sty, linewidth=linewidth, label='_nolegend_')
        else:
            plt.plot(x_values, points, line_sty, color=color, linewidth=linewidth, label='_nolegend_')
        txt_loc = np.where(x_values == np.min(x_values))
        plt.text(x_values[txt_loc], points[txt_loc], label, fontsize=8)

    def add_model_compare(self, root_name, color):
        frequencies, points = mdl.compare_nh3_model(root_name)
        self.labels_list_models.append((root_name + " : NH3 (new - old)"))
        plt.plot(self.f_or_w(frequencies), points, color=color)

    def add_all_other_data(self):
        for p in mdl.EXISTING_DATA:
            f, tb, tbe, c, m = p(self)
            plt.errorbar(self.f_or_w(f), tb, yerr=tbe,
                         fmt='.', color=c, marker=m)

    def f_or_w(self, frequencies):
        if self.wl:
            return (cst.c / np.array(frequencies)) * 1e-7
        else:
            return frequencies

    def show_prep(self, location='lower right'):
        plt.xscale('log', basex=2)
        if self.wl:
            tick_locations = [0.4, 0.6, 1, 2, 4, 6]
            plt.xticks(tick_locations, [loc for loc in tick_locations])
        else:
            tick_locations = np.arange(76) + 19
            plt.xticks(tick_locations, [str(loc) if loc % 10 == 0 else '' for loc in tick_locations])
        plt.legend(self.labels_list_models + self.labels_list, loc=location)

    def show(self):
        plt.show()
