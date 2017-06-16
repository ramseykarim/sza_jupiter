import numpy as np
from unpack import generate_names
import re


PATH = "/home/ramsey/Documents/Research/Jupiter/SZA/models/"
TRUE_MODEL_PATH = "/home/ramsey/Documents/Research/pyplanet/Output/"

WMAP_F = [22.44, 32.63, 40.50, 61.57, 92.71]
WMAP_TB = np.array([136.5, 149.0, 157.2, 168.0, 176.0])
WMAP_TBE = np.array([0.85, 0.77, 0.67, 0.59, 0.54]) * WMAP_TB / 100.

GIBSON_F = 28.5
GIBSON_TB = 142.9
GIBSON_TBE = 2.3


def generate_model(root_name):
    raw_data = np.loadtxt(PATH + root_name + "/" + root_name + ".dat")
    frequencies = raw_data[:, 0]
    points = raw_data[:, 1]
    return frequencies, points


def examine_model_text(name):
    matcher = re.compile("(\d+)d(\d+)_(\d+)pct")
    match = matcher.search(name)
    if match:
        ones_place = match.group(1)
        decimal_place = match.group(2)
        multiplier = match.group(3)
        p_cutoff = float(ones_place) + (float(decimal_place) / 10 ** len(decimal_place))
        multiplier = float(multiplier) / 100.
        return multiplier, p_cutoff
    else:
        raise NameError("The model signature was not found to be in the correct format.")


def examine_models(root_name_list):
    frequencies = 0
    temperature_list = []
    abundance_list = []
    cutoff_list = []
    matcher = re.compile("(\d+)d(\d+)_(\d+)pct")
    for name in root_name_list:
        match = matcher.search(name)
        if match:
            f, p = generate_model(name)
            temperature_list.append(p)
            frequencies = f
            ones_place = match.group(1)
            decimal_place = match.group(2)
            multiplier = match.group(3)
            p_cutoff = float(ones_place) + (float(decimal_place) / 10**len(decimal_place))
            cutoff_list.append(p_cutoff)
            multiplier = float(multiplier) / 100.
            abundance_list.append(multiplier)
        else:
            print "nope: " + name
    return frequencies, temperature_list, abundance_list, cutoff_list


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
    plotter.labels_list.append("KG 1976")
    return kg_f, kg_tb, kg_tbe, 'maroon', 'p'


def wmap(plotter):
    # wmap_f = [22.44, 32.63, 40.50, 61.57, 92.71]
    # wmap_tb = np.array([136.5, 149.0, 157.2, 168.0, 176.0])
    # wmap_tbe = np.array([0.85, 0.77, 0.67, 0.59, 0.54])
    # wmap_tbe = wmap_tb * wmap_tbe / 100.
    plotter.labels_list.append("WMAP 2011")
    return WMAP_F, WMAP_TB, WMAP_TBE, 'blue', '^'


def gibson(plotter):
    plotter.labels_list.append("Gibson 2005")
    return [28.5], [142.9], [2.3], 'red', 's'


def imke(plotter):
    f = [4.52, 5.49, 6.5, 7.5,
         8.5, 9.52, 10.46, 11.46,
         13.18, 14.21, 15.18, 16.21,
         17.38]
    tb = [247.5, 223.3, 207.6, 192.9,
          189.4, 178.9, 172.5, 167.3,
          167.1, 161.1, 156.8, 153.0,
          150.2]
    tbe = [7.4, 6.7, 6.2, 5.8,
           5.8, 5.6, 5.4, 5.3,
           5.3, 5.1, 4.9, 4.8,
           4.7]
    plotter.labels_list.append("de Pater 2016")
    return f, tb, tbe, 'orange', 'v'


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
                (model_frequencies, model_tb), deg=5, w=False, offset=False):
    # Fit to my data
    fit_space = np.where((model_frequencies > frequencies[frequencies.size - 1] - 1)
                         & (model_frequencies < frequencies[0] + 1))
    model_fit = np.polyfit(model_frequencies[fit_space], model_tb[fit_space], deg=deg)
    in_offset_range = True

    def model(x):
        output = np.zeros(x.size)
        for i, coefficient in enumerate(model_fit):
            output += coefficient * (x ** (deg - i))
        return output

    m_point = model(frequencies)
    mid_index = len(tb) / 2
    if not offset:
        error = tbe_slope
        offset = m_point[mid_index] - tb[mid_index]
    else:
        offset = 0.
        error = tbe_offset
    # if np.abs(offset) > np.min(tbe_offset):
    #     in_offset_range = False
    chi_squared_value = m_point - tb - offset
    df = len(chi_squared_value)  # Degrees of Freedom (number of points)
    chi_squared_value = np.sum((chi_squared_value ** 2.) / (error ** 2.))  # TODO Change to tbe_slope/offset to see diff
    # Fit to WMAP:
    if w:
        near_wmap = np.concatenate([model_tb[0:1], model_tb[-3:]])
        df += len(near_wmap)
        wmap_ch_sq = np.sum(((near_wmap - WMAP_TB) ** 2.) / WMAP_TBE**2.)
        chi_squared_value += wmap_ch_sq

        # Fitting Gibson
        f_loc = np.where(model_frequencies == GIBSON_F)
        df += 1
        gibson_ch_sq = ((GIBSON_TB - model_tb[f_loc]) / GIBSON_TBE)**2.
        chi_squared_value += gibson_ch_sq

    df -= (1 + (1 - offset))

    return chi_squared_value/df, in_offset_range


def create_parameter_space((s1, e1, step1), (s2, e2, step2), info_tuple, model_path, w=False):
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
        chi_squared_value, in_range = chi_squared(info_tuple, model_iterator.next(), w=w)
        chi_squared_values.append((chi_squared_value, in_range))
    return params, chi_squared_values


def create_parameter_space_1d(param_kws, info_tuple, model_path, w=False, offset=None):
    if type(param_kws) is tuple:
        params = np.arange(param_kws[0], param_kws[1], param_kws[2])
    else:
        params = param_kws
    model_path = TRUE_MODEL_PATH + model_path
    model_iterator = ModelGrid(model_path)
    chi_squared_values = []
    while model_iterator.has_next():
        chi_squared_value, in_range = chi_squared(info_tuple, model_iterator.next(), w=w, offset=offset)
        chi_squared_values.append(chi_squared_value)
    return params, np.array(chi_squared_values, dtype=float)


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
