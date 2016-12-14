import numpy as np
import channel as ch
import random
import meanmethods as mm


def shuffle(array, error):
    order = range(len(array))
    random.shuffle(order)
    array_c = np.array([])
    error_c = np.array([])
    for i in order:
        array_c = np.append(array_c, array[i])
        error_c = np.append(error_c, error[i])
    return array_c, error_c


def progressive_means(array, error, n):
    method_results = mm.MethodContainer(array.size - 2)
    np.seterr(all='ignore')
    for i in range(n):
        a, e = shuffle(array, error)
        for j in range(2, a.size):
            methods = mm.AverageMethods(a[:j], e[:j])
            for k, f in enumerate(methods.functions):
                method_results.results[k][j - 2] += f()
    for i in range(len(method_results.results)):
        method_results.results[i] /= n
    return method_results


def bootstrap_flux_average(array, error, frequency, n):
    half = len(array) / 2
    results = mm.BootstrapAverageContainer(frequency)
    for i in range(n):
        a, e = shuffle(array, error)
        average_container = mm.AverageMethods(a[:half], e[:half])
        results.append(average_container.method_5_a())
    return results


class Stats:
    def __init__(self, channel_list):
        assert isinstance(channel_list, list) and isinstance(channel_list[0], ch.Channel)
        self.f_e_f_triples = [(np.copy(c.flux), np.copy(c.error), c.frequency_ghz) for c in channel_list]
        self.results_error = []
        self.results_averages = []

    def execute_error_testing(self, n):
        for f, e, v in self.f_e_f_triples:
            print '.',
            self.results_error.append(progressive_means(f, e, n))
        print '!'

    def execute_flux_testing(self, n):
        for f, e, v in self.f_e_f_triples:
            print '.',
            self.results_averages.append(bootstrap_flux_average(f, e, v, n))
        print '!'
        return self.results_averages


