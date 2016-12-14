import numpy as np
import channel as ch
import unpack as up
import plotting as pltt
import matplotlib.pyplot as plt
import sync_cmb as sc


class FlexChannel(ch.Channel):
    def __init__(self, rectangle, frequency_ghz):
        ch.Channel.__init__(self, rectangle, frequency_ghz)

    def average_choose(self, number):
        if number == 1:
            self.average()
        elif number == 2:
            self.average_error_jy()
        elif number == 3:
            self.average_n_hat()
        elif number == 4:
            self.average_sigma_hat()
        else:
            raise ValueError(str(number) + " was not an option")

    def average_error_jy(self):
        """
        This is the same as "new" but with error assumed in jy
        :return: The averages
        """
        assert isinstance(self.error, np.ndarray)
        assert isinstance(self.flux, np.ndarray)
        weights = 1. / (self.error * self.flux) ** 2.
        variance_sq = 1. / np.nansum(weights)
        self.flux_avg = np.nansum(weights * self.flux) * variance_sq
        acc_error_sq = np.nansum(weights * ((self.flux - self.flux_avg) ** 2.)) * variance_sq
        total_error_sq_inv = (1. / acc_error_sq) + (1. / variance_sq)
        self.error_avg = 1. / np.sqrt(total_error_sq_inv)

    def average_n_hat(self):
        """
        Only uses n hat as error
        :return: The averages
        """
        assert isinstance(self.error, np.ndarray)
        assert isinstance(self.flux, np.ndarray)
        weights = 1. / self.error ** 2.
        variance_sq = 1. / np.nansum(weights)
        self.flux_avg = np.nansum(weights * self.flux) * variance_sq
        self.error_avg = np.sqrt(variance_sq)

    def average_sigma_hat(self):
        """
        Only uses sigma hat as error
        :return: The averages
        """
        assert isinstance(self.error, np.ndarray)
        assert isinstance(self.flux, np.ndarray)
        weights = 1. / self.error ** 2.
        variance_sq = 1. / np.nansum(weights)
        self.flux_avg = np.nansum(weights * self.flux) * variance_sq
        acc_error_sq = np.nansum(weights * (self.flux - self.flux_avg) ** 2.) * variance_sq
        self.error_avg = np.sqrt(acc_error_sq)


def flex_unpack_copss(csv_file, frequency_list_ghz):
    raw_all = np.genfromtxt(csv_file)
    date_array = raw_all[:, 0]
    channel_list = []
    for i in range(15):
        channel_list.append(FlexChannel(raw_all[:, i * 2 + 1:i * 2 + 3], frequency_list_ghz[i]))
    return date_array, channel_list


class FlexUnpack(up.Unpack):
    def __init__(self):
        up.Unpack.__init__(self)
        self.dates_copss_array, self.channel_obj_list = \
            flex_unpack_copss(up.PATH + "TXTs/jupiter_fluxes_flagged_USETHIS.txt",
                              self.frequency_list_ghz)
        self.number = -1

    def assign(self, number):
        self.number = number
        return self

    def average_adj(self):
        for channel in self.channel_obj_list:
            channel.average_choose(self.number)


def plot_mean_bars(mean, frequency, color='k'):
    plt.plot([frequency - 0.03, frequency + 0.03], [mean, mean], color=color)


def add_best_mean_bars(unpacker, result_list):
    pairs = [(channel.tb, channel.frequency_ghz) for channel in unpacker.channel_obj_list]
    [plot_mean_bars(m, f, color='blue') for m, f in pairs]
    [plot_mean_bars(sc.synchrotron_cmb(c.mean(), c.x[0]), c.x[0]) for c in result_list]


class FlexPlotting(pltt.Plotting):
    def __init__(self, wl=False):
        pltt.Plotting.__init__(self, wl=wl)
        plt.autoscale()

    def add_my_data(self, unpacker):
        self.custom_label(unpacker.number)
        tb = [channel.tb for channel in unpacker.channel_obj_list]
        tb_err = [channel.tb_error for channel in unpacker.channel_obj_list]
        x_axis = np.array(unpacker.frequency_list_ghz) + (float(unpacker.number) - 2.5) * 0.01
        plt.errorbar(x_axis, tb, yerr=tb_err, fmt='.')

    def add_flux_scatter(self, result_list):
        x = np.array([])
        y = np.array([])
        t = np.array([])
        means = np.array([])
        errors = np.array([])
        for container in result_list:
            x = np.append(x, np.array(container.x))
            y = np.append(y, np.array(container.y))
            t = np.append(t, np.array(container.t))
            means = np.append(means, container.mean())
            errors = np.append(errors, container.s_d_o_m())
        fs = np.array([f.x[0] for f in result_list])
#        plt.plot(x, y, ',', color='k')
#        plt.errorbar((fs + 0.002), means, yerr=errors, fmt='.', color='k')
        hi_error = sc.synchrotron_cmb_units_np(means + errors, fs)
        lo_error = sc.synchrotron_cmb_units_np(means - errors, fs)
        means = sc.synchrotron_cmb_units_np(means, fs)
        hi_error, lo_error = np.abs(hi_error - means), np.abs(lo_error - means)
        error_bounds = np.array([hi_error, lo_error])
        errors = np.max(error_bounds, axis=0)
        plt.plot(x, t, ',', color='k')
        self.labels_list.append("Bootstrapping")
        plt.errorbar((fs + 0.002), means, yerr=errors, fmt='.', color='k')
        self.labels_list.append("Bootstrap mean")

    def show_prep(self):
        plt.legend(self.labels_list_models + self.labels_list, loc='lower right')

    def show(self):
        plt.show()

    def custom_label(self, number):
        if number == 1:
            self.labels_list.append("$\sigma_{T}$ w/ jy")
        elif number == 2:
            self.labels_list.append("$\sigma_{T}$ not jy")
        elif number == 3:
            self.labels_list.append("Therm only")
        elif number == 4:
            self.labels_list.append("Ensemble")
        else:
            raise ValueError(str(number) + " was not an option")

