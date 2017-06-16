import numpy as np
import sys
import unpack as up
import channel as ch
import random
import matplotlib.pyplot as plt

FREQUENCIES = [34.688000, 34.188000, 33.688000,
               33.188000, 32.688000, 32.188000,
               31.688000, 31.188000, 30.688000,
               30.188000, 29.688000, 29.188000,
               28.688000, 28.188000, 27.688000]


def even_odd(u, c):
    array = c.flux
    return np.where(np.arange(len(array)) % 2 == 0)


def first_last(u, c):
    array = c.flux
    quarter_length = len(array) / 4
    index_array = np.arange(len(array))
    return np.where((index_array > 3 * quarter_length) | (index_array < quarter_length))


def first_half(u, c):
    array = c.flux
    index_array = np.arange(len(array))
    half_length = len(array) / 2
    return_array = np.where(index_array < half_length)
    return return_array


def j_elevation(u, c):
    array = u.get_j_el()
    return np.where(array < np.mean(array) - 0.5 * np.std(array))


def m_elevation(u, c):
    array = u.get_m_el()
    return np.where(array < np.median(array) - np.std(array))


def debug(u, c):
    array = c.flux
    return np.arange(len(array))


JACKKNIFE_TESTS = {
    "Even-Odd": (even_odd, 'o', "blue"),
    "First-Last": (first_last, '^', "green"),
    "Jupiter El": (j_elevation, 's', "red"),
    # "Mars El": (m_elevation, 'D', "maroon"),
    # "Half": (first_half, 'p', "cyan"),
    "Debug": (debug, '.', "orange"),
}


class JackknifeUnpack(up.Unpack):
    def __init__(self):
        up.Unpack.__init__(self)
        self.channel_type = JackknifeChannel
        self.test = lambda x: sys.stdout.write("You didn't assign a test!")
        up.Unpack.prepare(self)
        self.original_channel_objects = [c.copy() for c in self.channel_obj_list]

    def accept_test(self, test):
        self.prepare()
        self.test = test
        return self

    def prelim_adjust(self):
        for channel in self.channel_obj_list:
            channel.pick_data(self)
        return self

    def prepare(self):
        self.channel_obj_list = [c.copy() for c in self.original_channel_objects]
        return self


class JackknifeChannel(ch.Channel):
    def __init__(self, rectangle, frequency_ghz):
        ch.Channel.__init__(self, rectangle, frequency_ghz)

    def pick_data(self, unpacker):
        indices = unpacker.test(unpacker, self)
        self.flux, self.error = self.flux[indices], self.error[indices]


class BootstrapAverageContainer:
    def __init__(self, frequencies, unpacker):
        self.u = unpacker
        self.frequencies = frequencies
        self.f = np.array([])
        self.t = np.array([])

    def cycle(self, offset):
        self.u.accept_test(shuffle)
        for i in range(500):
            if i % 5 == 0:
                msg = "Bootstrapping {0} % Complete".format(i/5 + 1).strip()
                sys.stdout.write(msg)
                sys.stdout.write("\r")
                sys.stdout.flush()
            self.u.prepare().adjust()
            tb, unused, also_unused = self.u.get_temperatures()
            self.append(tb - offset)

    def append(self, temperatures):
        self.f = np.concatenate([self.f, self.frequencies])
        self.t = np.concatenate([self.t, temperatures])


def shuffle(u, c):
    array = c.flux
    order = range(len(array))
    random.shuffle(order)
    return order[:len(order)/2]

