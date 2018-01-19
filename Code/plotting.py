import matplotlib.pyplot as plt
import models as mdl
import unpack as up
import datacontainer as dc
import scipy.constants as cst
import numpy as np


class Plotting:
    def __init__(self, wl=False):
        self.wl = wl

    def prep_data_plot(self):
        plt.xlabel("$\nu$ (GHz)" if not self.wl else "$\lambda$ (cm)")
        plt.ylabel("$T_{b}$ (K)")
        plt.title("$T_{b}$ vs. " + "$\\nu$" if not self.wl else "$\lambda$")
        plt.xlim(20, 95)
        plt.ylim(120, 220)

    def add_my_data(self, unpacker, ax=None, psize=1.2):
        if isinstance(unpacker, dc.DataContainer):
            x_axis, tb, tb_rel_err, tb_ens_err = unpacker.yield_model_comp_info()
        else:
            tb = [channel.tb for channel in unpacker.channel_obj_list]
            tb_rel_err = np.array([channel.tb_error_slope for channel in unpacker.channel_obj_list])
            tb_ens_err = np.array([channel.tb_error_offset for channel in unpacker.channel_obj_list])
            x_axis = unpacker.frequency_list_ghz if not self.wl else [c.wavelength_cm for c in
                                                                      unpacker.channel_obj_list]
        combined_error = np.sqrt(tb_rel_err ** 2. + tb_ens_err ** 2.)
        if ax is not None:
            thing = ax
        else:
            thing = plt
        thing.errorbar(x_axis, tb, yerr=tb_ens_err, fmt='o', markersize=5.,
                       color='grey', mec='black', mfc='white', mew=psize, label='_nolegend_')
        thing.errorbar(x_axis, tb, yerr=tb_rel_err, fmt='o', capsize=4., markersize=5.,
                       color='black', mec='black', mfc='white', mew=psize, label='CARMA')

    def add_model_plot(self, model_path, color='k', line_sty='-', label=None, linewidth=2.0, ax=None):
        if label is None:
            label = "_nolegend_"
        frequencies, points = mdl.generate_model(model_path)
        x_values = self.f_or_w(frequencies)
        if ax is not None:
            thing = ax
        else:
            thing = plt
        if color.lower() == "default":
            thing.plot(x_values, points, line_sty, linewidth=linewidth, label=label)
        else:
            thing.plot(x_values, points, line_sty, color=color, linewidth=linewidth, label=label)
            # txt_loc = np.where(x_values == np.min(x_values))
            # plt.text(x_values[txt_loc], points[txt_loc], label, fontsize=8)

    def add_all_other_data(self, ax=None, psize=3.):
        if ax is not None:
            thing = ax
        else:
            thing = plt
        for p in mdl.EXISTING_DATA:
            f, tb, tbe, c, m, l = p()
            thing.errorbar(self.f_or_w(f), tb, yerr=tbe, capsize=3.,
                           fmt='.', color=c, marker=m, label=l,
                           elinewidth=.5, markersize=psize)

    def f_or_w(self, frequencies):
        if self.wl:
            return (cst.c / np.array(frequencies)) * 1e-7
        else:
            return frequencies

    def show_prep(self, location='lower right', ncol=1, legend=True):
        plt.xscale('log', basex=2)
        if self.wl:
            tick_locations = [0.4, 0.6, 1, 2, 4, 6]
            plt.xticks(tick_locations, [loc for loc in tick_locations])
        else:
            tick_locations = np.arange(76) + 19
            plt.xticks(tick_locations, [str(loc) if loc % 10 == 0 else '' for loc in tick_locations])
        if legend:
            plt.legend(loc=location, ncol=ncol)

    def show(self):
        plt.show()
