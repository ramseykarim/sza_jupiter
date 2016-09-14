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


class Stats:
    def __init__(self, channel_list):
        assert isinstance(channel_list, list) and isinstance(channel_list[0], ch.Channel)
        self.f_e_pairs = [(np.copy(c.flux), np.copy(c.error)) for c in channel_list]
        self.results = []

    def execute_testing(self, n):
        for f, e in self.f_e_pairs:
            print '.',
            self.results.append(progressive_means(f, e, n))
        print '!'


