import matplotlib.pyplot as plt
import models as mdl
import scipy.constants as cst


class Plotting:
    def __init__(self, wl=False):
        self.wl = wl
        self.labels_list = []
        self.labels_list_models = []
        plt.xlabel("Frequency (GHz)" if not self.wl else "Wavelength (cm)")
        plt.ylabel("$T_{b}$ (K)")
        plt.title("$T_{b}$, synchrotron & CMB corrected")

    def add_my_data(self, unpacker):
        self.labels_list.append("SZA")
        tb = [channel.tb for channel in unpacker.channel_obj_list]
        tb_err = [channel.tb_error for channel in unpacker.channel_obj_list]
        x_axis = unpacker.frequency_list_ghz if not self.wl else [c.wavelength_cm for c in unpacker.channel_obj_list]
        plt.errorbar(x_axis, tb, yerr=tb_err, fmt='o', color="red")

    def add_model_plot(self, root_name, color):
        frequencies, points = mdl.generate_model(root_name)
        self.labels_list_models.append(root_name)
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

    def show(self):
        plt.legend(self.labels_list_models + self.labels_list)
        plt.show()
