import numpy as np

"""
All methods are designed to take two arrays:
    :param array
        is to be the flux array, in units of janskys
    :param error
        is to be the error array, as a fraction of the flux

All return either the averaged flux in janskys or the averaged error ALSO in janskys
"""


def coerce_to_np(array):
    if not isinstance(array, np.ndarray):
        return np.array(array)
    else:
        return array


class AverageMethods:
    def __init__(self, array, error):
        self.a = coerce_to_np(array)
        self.e = coerce_to_np(error)
        self.weights = 1. / (self.a * self.e)
        self.w_sum = np.nansum(self.weights)
        self.f_avg = self.flux_avg()
        self.f_avg_2 = 0.
        self.functions = [self.method_1_a, self.method_2_a,
                          self.method_1_e, self.method_2_e,
                          self.method_3_e, self.method_4_e,
                          self.method_5_e]

    def flux_avg(self):
        """
        :return: Normal weighted average flux
        """
        weighted_flux = 1. / self.e
        sum_flux = np.nansum(weighted_flux)
        avg_flux = sum_flux / self.w_sum
        return avg_flux

    def method_1_a(self):
        """
        :return: Std Deviation
        """
        array = self.a
        return np.std(array)

    def method_2_a(self):
        """
        :return: Normal weighted flux average
        """
        return self.f_avg

    def method_1_e(self):
        """
        :return: Systematic-ignorant thermal noise reduction
        """
        initial = np.sqrt(np.nansum(self.weights ** 2.))
        return 1. / initial

    def method_2_e(self):
        """
        :return: Weighted standard deviation, essentially -- from poster
        """
        avg_flux = self.f_avg
        initial = np.nansum(self.weights * (self.a - avg_flux) ** 2.) / self.w_sum
        return np.sqrt(initial)

    def method_3_e(self):
        """
        :return: What we wrote down in a meeting
        """
        avg_flux = self.f_avg
        initial = np.sqrt(np.nansum(self.weights ** 2. * (self.a - avg_flux) ** 2.))
        return 1. / initial

    def method_4_e(self):
        """
        :return: Some bad idea I had
        """
        avg_flux = self.f_avg
        initial = np.sqrt(np.nansum((self.weights / (self.a - avg_flux)) ** 2.))
        return 1. / initial

    def method_5_a(self):
        """
        :return: Different flux averaging technique, using squared weight
        """
        weights2 = self.weights ** 2.
        initial = np.nansum(weights2 * self.a) / np.nansum(weights2)
        self.f_avg_2 = initial
        return initial

    def method_5_e(self):
        """
        :return: Best error averaging. But rather small error?
        """
        avg_flux = self.f_avg_2
        weights2 = self.weights ** 2.
        variance_sq = 1. / np.nansum(weights2)
        acc_error_sq = np.nansum(weights2 * (self.a - avg_flux)**2.) * variance_sq
        total_error_sq_inv = (1. / acc_error_sq) + (1. / variance_sq)
        return 1. / np.sqrt(total_error_sq_inv)


class MethodContainer:
    def __init__(self, size):
        self.method_1_a = np.zeros(size)
        self.method_1_e = np.zeros(size)
        self.method_2_a = np.zeros(size)
        self.method_2_e = np.zeros(size)
        self.method_3_e = np.zeros(size)
        self.method_4_e = np.zeros(size)
        self.method_5_e = np.zeros(size)
        self.results = [self.method_1_a, self.method_2_a,
                        self.method_1_e, self.method_2_e,
                        self.method_3_e, self.method_4_e,
                        self.method_5_e]
