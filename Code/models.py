import numpy as np
from unpack import generate_names


PATH = "/home/ramsey/Documents/Research/Jupiter/SZA/models/"
TRUE_MODEL_PATH = "/home/ramsey/Documents/Research/pyplanet/Output/"


def generate_model(root_name):
    raw_data = np.loadtxt(PATH + root_name + "/" + root_name + ".dat")
    frequencies = raw_data[:, 0]
    points = raw_data[:, 1]
    return frequencies, points


def compare_nh3_model(root_name):
    new_s = "_newnh3"
    raw_data_old = np.loadtxt(PATH + root_name + "/" + root_name + ".dat")
    frequencies = raw_data_old[:, 0]
    points_old = raw_data_old[:, 1]
    raw_data_new = np.loadtxt(PATH + root_name + new_s + "/" + root_name + new_s + ".dat")
    points_new = raw_data_new[:, 1]
    return frequencies, (points_new - points_old)


def kg(plotter):
    kg_f = [20.050, 20.295, 20.735, 20.970, 21.225, 21.550, 22.055, 22.220, 22.370, 22.945, 23.410, 23.840, 24.100]
    kg_tb = [144.9, 145.6, 144.0, 143.0, 146.0, 142.2, 142.1, 137.6, 137.1, 138.6, 139.1, 141.6, 140.8]
    kg_tbe = [5.8, 5.8, 5.8, 5.7, 5.8, 5.7, 5.7, 5.5, 5.5, 5.5, 5.6, 5.7, 5.6]
    plotter.labels_list.append("KG")
    return kg_f, kg_tb, kg_tbe


def wmap(plotter):
    wmap_f = [22.44, 32.63, 40.50, 61.57, 92.71]
    wmap_tb = np.array([136.5, 149.0, 157.2, 168.0, 176.0])
    wmap_tbe = np.array([0.85, 0.77, 0.67, 0.59, 0.54])
    wmap_tbe = wmap_tb * wmap_tbe / 100.
    plotter.labels_list.append("WMAP")
    return wmap_f, wmap_tb, wmap_tbe


def gibson(plotter):
    plotter.labels_list.append("Gibson")
    return [28.5], [142.9], [2.3]


def imke(plotter):
    f = [13.18, 14.21, 15.18, 16.21, 17.38, 8.5, 9.52, 10.46, 11.46, 4.52, 5.49, 6.5, 7.5]
    tb = [159.5, 154.5, 150.6, 146.4, 143.5, 187.7, 177.5, 172.4, 165.8, 247.5, 223.3, 207.6, 192.9]
    tbe = [4.8, 4.6, 4.5, 4.4, 4.3, 5.7, 5.4, 5.2, 5., 7.4, 6.7, 6.2, 5.8]
    plotter.labels_list.append("de Pater")
    return f, tb, tbe


EXISTING_DATA = [kg, wmap, gibson, imke]


def old_chi_sq((frequencies, tb, tbe_small, tbe_large, model_name)):
    m_freq, m_point = generate_model(model_name)

    deg = 3
    model_fit = np.polyfit(m_freq, m_point, deg)

    def model(x):
        output = np.zeros(x.shape)
        for i, coefficient in enumerate(model_fit):
            output += coefficient * (x ** (deg - i))
        return output

    m_point = model(frequencies)
    mid_index = len(tb) / 2
    offset = m_point[mid_index] - tb[mid_index]
    if np.abs(offset) > np.min(tbe_large):
        return False
    chi_squared_value = m_point - tb - offset
    chi_squared_value = np.sum((chi_squared_value ** 2.) / tbe_small)
    return chi_squared_value


def chi_squared((frequencies, tb, tbe_slope, tbe_offset),
                (model_frequencies, model_tb), deg=5):
    fit_space = np.where((model_frequencies > frequencies[frequencies.size - 1]) & (model_frequencies < frequencies[0]))
    model_fit = np.polyfit(model_frequencies[fit_space], model_tb[fit_space], deg=deg)
    in_offset_range = True

    def model(x):
        output = np.zeros(x.size)
        for i, coefficient in enumerate(model_fit):
            output += coefficient * (x ** (deg - i))
        return output

    m_point = model(frequencies)
    mid_index = len(tb) / 2
    offset = m_point[mid_index] - tb[mid_index]
    if np.abs(offset) > np.min(tbe_offset):
        in_offset_range = False
    chi_squared_value = m_point - tb - offset
    chi_squared_value = np.sum((chi_squared_value ** 2.) / tbe_slope)
    return chi_squared_value, in_offset_range


def create_parameter_space((s1, e1, step1), (s2, e2, step2), info_tuple, model_path):
    param1 = np.arange(s1, e1, step1)
    param2 = np.arange(s2, e2, step2)
    params = []
    for p1 in param1:
        for p2 in param2:
            params.append((p1, p2))
    model_path = TRUE_MODEL_PATH + model_path
    model_iterator = ModelGrid(model_path)
    chi_squared_values = []
    while model_iterator.has_next():
        chi_squared_value, in_range = chi_squared(info_tuple, model_iterator.next())
        chi_squared_values.append((chi_squared_value, in_range))
    return params, chi_squared_values


class ModelGrid:
    def __init__(self, path):
        self.path = path
        self.file_names = generate_names(path)
        self._index = -1

    def next(self):
        if not self.has_next():
            raise IndexError
        self._index += 1
        file_name = self.path + self.file_names[self._index]
        data = np.loadtxt(file_name)
        return data[:, 0], data[:, 1]

    def has_next(self):
        return self._index < len(self.file_names) - 1

    def current(self):
        self._index -= 1
        return self.next()

    def reinitialize(self):
        self._index = 1

    def get_index(self):
        return self._index
