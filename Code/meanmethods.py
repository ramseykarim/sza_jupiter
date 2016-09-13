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
        self.last = 0
        self.functions = [self.method_1_a, self.method_1_e,
                          self.method_2_a, self.method_2_e,
                          self.method_3_a, self.method_3_e,
                          self.method_4_a, self.method_4_e]

    def method_1_a(self):
        array = self.a
        return np.std(array)

    def method_1_e(self):
        initial = np.sqrt(np.nansum(self.weights ** 2.))
        return 1. / initial

    def method_2_a(self):
        weighted_flux = 1. / self.e
        sum_flux = np.nansum(weighted_flux)
        avg_flux = sum_flux / self.w_sum
        self.last = avg_flux
        return avg_flux

    def method_2_e(self):
        avg_flux = self.last
        initial = np.nansum(self.weights * (self.a - avg_flux) ** 2.) / self.w_sum
        return np.sqrt(initial)

    def method_3_a(self):
        self.last = self.method_2_a()
        return self.last

    def method_3_e(self):
        avg_flux = self.last
        initial = np.sqrt(np.nansum(self.weights ** 2. * (self.a - avg_flux) ** 2.))
        return 1. / initial

    def method_4_a(self):
        self.last = self.method_2_a()
        return self.last

    def method_4_e(self):
        avg_flux = self.last
        initial = np.sqrt(np.nansum((self.weights * (self.a - avg_flux) / avg_flux) ** 2.))
        return 1. / initial


class MethodContainer:
    def __init__(self, size):
        self.method_1_a = np.zeros(size)
        self.method_1_e = np.zeros(size)
        self.method_2_a = np.zeros(size)
        self.method_2_e = np.zeros(size)
        self.method_3_a = np.zeros(size)
        self.method_3_e = np.zeros(size)
        self.method_4_a = np.zeros(size)
        self.method_4_e = np.zeros(size)
        self.results = [self.method_1_a, self.method_1_e,
                        self.method_2_a, self.method_2_e,
                        self.method_3_a, self.method_3_e,
                        self.method_4_a, self.method_4_e]
