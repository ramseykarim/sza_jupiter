import numpy as np
from unpack import generate_names
from numbers import Number
import re

PATH = "/home/ramsey/Documents/Research/Jupiter/SZA/models/"
TRUE_MODEL_PATH = "/home/ramsey/Documents/Research/pyplanet/Output/"

WMAP_F = np.array([22.44, 32.63, 40.50, 61.57, 92.71])
WMAP_TB = np.array([136.5, 149.0, 157.2, 168.0, 176.0])
WMAP_TBE = np.array([0.85, 0.77, 0.67, 0.59, 0.54]) * WMAP_TB / 100.

GIBSON_F = 28.5
GIBSON_TB = 142.9
GIBSON_TBE = 2.3


def physical_model(p, original_nh3, top, humidity, deep):
    """
    Assumes p, nh3 sorted by p
    :param p:
    :param original_nh3:
    :param top:
    :param humidity:
    :param deep:
    :return:
    """
    return


def examine_model(model_name):
    matcher = re.compile("top(\d{3})_rh(\d{3})_deep(\d{3})\.dat")
    match = matcher.search(model_name)
    top = int(match.group(1))
    rh = int(match.group(2))
    deep = int(match.group(3))
    return top, rh, deep


def find_saturation_curve(p, original_nh3):
    """
    The plan is to never change the order of the layers.
    :param p: p
    :param original_nh3: nh3
    :return: Indices of the cloud top and bottom, as well as the
        fit parameters of the saturation curve.
    """
    sorted_order = np.argsort(p)
    i = 0
    while original_nh3[sorted_order[i]] == original_nh3[sorted_order[i + 1]]:
        i += 1
    cloud_top = sorted_order[i]
    while original_nh3[sorted_order[i]] < original_nh3[sorted_order[i + 1]]:
        i += 1
    cloud_bottom = sorted_order[i]
    lo_lim = min(cloud_top, cloud_bottom)
    hi_lim = max(cloud_top, cloud_bottom) + 1
    p_curve = p[lo_lim:hi_lim]
    nh3_curve = original_nh3[lo_lim:hi_lim]
    deg = 8
    curve_fit = np.polyfit(p_curve, np.log(nh3_curve), deg=deg)
    return cloud_top, cloud_bottom, curve_fit


def physical_tcm_model_prototype(top, rh, deep):
    """
    The plan is to never change the order of the layers!
    The algorithm is as follows:
        1) Fit the saturation curve. This yields the original
            NH3 cloud top and bottom. The old cloud bottom
            pressure is P_sat.
        2) Save the modified TOP and DEEP (~1 bar) abundances.
        3) Apply DEEP to 8 bar < P < 1.7 bar. 1.7 bar is where
            the previous cloud's depletion ends.
        4) Step up from 1.7 bar, applying the new "DEEP" ~1bar
            abundance (1.7 bar works) until (P/RH = P_eff) < P_sat.
        5) Above this new cloud base, apply the shifted fit,
            keeping in mind that the shift is logarithmic (P/RH)
            and that an exponential must be applied to the fit.
            Do this until the next NH3 value would fall below
            the new TOP abundance.
        6) After this new cloud top, apply the TOP abundance.
    Note that P > 8 bar is left alone, per Galileo results.
    :param top: Multiplier for the uppermost, constant dry section.
        Applied at P < SaturationCurve.
    :param rh: Multiplier for relative humidity. This is applied by
        shifting the saturation curve up or down, changing the
        cloud base location.
    :param deep: Multiplier for the deep abundance. Applied between
        8 bar > P > SaturationCurve.
    :return:
    """
    base_model = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulvla72',
                               skip_header=1)
    p = base_model[:, 2]
    nh3 = base_model[:, 6]
    cloud_top, cloud_bottom, curve_fit = find_saturation_curve(p, nh3)
    galileo_lim = 8
    nh4sh_lim = 1.7
    nh3[np.where((galileo_lim > p) & (p >= nh4sh_lim))] *= deep
    upper = np.where(p < nh4sh_lim)
    p_upper = p[upper]
    upper_i_sorted = upper[0][np.argsort(np.reciprocal(p_upper))]
    i = 0
    top_abundance = float(nh3[np.where(p == np.min(p))]) * top
    nh3[upper_i_sorted[i]] *= deep
    deep_abundance = nh3[upper_i_sorted[i]]
    p_sat_0 = p[cloud_bottom]
    # Used to use p[upper_i_sorted[i]]/rh > p_sat_0
    while p[upper_i_sorted[i]] > p_sat_0 and nh3[upper_i_sorted[i + 1]] >= nh3[upper_i_sorted[i]]:
        nh3[upper_i_sorted[i + 1]] = deep_abundance
        i += 1
    within_cloud = True
    new_cloud_bottom = i
    while within_cloud:
        # Used to use p[upper_i_sorted[i]]/rh
        updated_nh3 = np.exp(polynomial_single(curve_fit, p[upper_i_sorted[i]])) * rh
        if updated_nh3 > top_abundance:
            if updated_nh3 > deep_abundance:
                nh3[upper_i_sorted[i]] = deep_abundance
            else:
                nh3[upper_i_sorted[i]] = updated_nh3
            i += 1
        else:
            if i == new_cloud_bottom:
                print "<<================================>>"
                print "Something went wrong. Diagnostics:"
                print "TOP", top
                print "RH", rh
                print "DEEP", deep
                print "i", i
                print "current pressure", p[upper_i_sorted[i]]
                print "top abundance", top_abundance
                print "deep abundance", deep_abundance
                print "prev abundance", nh3[upper_i_sorted[i - 1]]
                print "attempted abundance", updated_nh3
                print "<<================================>>"
                raise IndexError("Your ammonia ice cloud has dropped below another cloud.")
            within_cloud = False
    nh3[upper_i_sorted[i:]] = top_abundance
    return p, nh3


def default_model():
    base_model = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulvla72',
                               skip_header=1)
    p = base_model[:, 2]
    nh3 = base_model[:, 6]
    return p, nh3


def generate_model(model_path):
    raw_data = np.loadtxt(model_path)
    frequencies = raw_data[:, 0]
    points = raw_data[:, 1]
    return frequencies, points


def kg():
    kg_f = [20.050, 20.295, 20.735, 20.970, 21.225, 21.550, 22.055, 22.220, 22.370, 22.945, 23.410, 23.840, 24.100]
    kg_tb = [144.9, 145.6, 144.0, 143.0, 146.0, 142.2, 142.1, 137.6, 137.1, 138.6, 139.1, 141.6, 140.8]
    kg_tbe = [5.8, 5.8, 5.8, 5.7, 5.8, 5.7, 5.7, 5.5, 5.5, 5.5, 5.6, 5.7, 5.6]
    return kg_f, kg_tb, kg_tbe, 'maroon', 'p', "KG 1976"


def wmap():
    return WMAP_F, WMAP_TB, WMAP_TBE, 'blue', '^', "WMAP 2011"


def gibson():
    # return [28.5], [142.9], [2.3], 'red', 's', "Gibson 2005"
    return [GIBSON_F], [GIBSON_TB], [GIBSON_TBE], 'red', 's', "Gibson 2005"


def imke():
    f = [4.52, 5.49, 6.5, 7.5,
         8.5, 9.52, 10.46, 11.46,
         13.18, 14.21, 15.18, 16.21,
         17.38
         ]

    tb = [247.5, 223.3, 207.6, 192.9,
          183.3, 177.7, 172.8, 167.6,
          164.0, 158.7, 155.0, 151.3,
          148.5
          ]

    tbe = [7.4, 6.7, 6.2, 5.8,
           5.499, 5.331, 5.184, 5.028,
           4.92, 4.761, 4.65, 4.539,
           4.455
           ]
    return f, tb, tbe, 'orange', 'v', "de Pater 2016"


def imke_old():
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
    return f, tb, tbe, 'orange', 'v', "de Pater 2016"


EXISTING_DATA = [kg, wmap, gibson, imke]


def chi_squared((frequencies, tb, tbe_slope, tbe_offset),
                (model_frequencies, model_tb), deg=5, all_others=False):
    # Fit to my data
    fit_space = np.where((model_frequencies > frequencies[frequencies.size - 1] - 1)
                         & (model_frequencies < frequencies[0] + 1))
    model_fit = np.polyfit(model_frequencies[fit_space], model_tb[fit_space], deg=deg)

    def model(x):
        return polynomial(model_fit, x)

    bad_freq = 33.188
    freqs_limited_33d188 = frequencies[np.where(frequencies != bad_freq)]
    matching_points = model(freqs_limited_33d188)
    tb_limited_33d188 = tb[np.where(frequencies != bad_freq)]
    mid_index = len(tb) / 2
    if not all_others:
        # Slope only
        error = tbe_slope[np.where(frequencies != bad_freq)]
        offset = matching_points[mid_index] - tb_limited_33d188[mid_index]
    else:
        # Slope and offset
        error = tbe_offset[np.where(frequencies != bad_freq)]
        offset = 0.
    chi_squared_value = matching_points - tb_limited_33d188 - offset
    df = len(chi_squared_value)
    chi_squared_value = np.sum((chi_squared_value ** 2.) / (error ** 2.))
    if all_others:
        # Fit to Gibson
        gibson_f_loc = np.where(model_frequencies == GIBSON_F)
        df += 1
        gibson_ch_sq = ((GIBSON_TB - model_tb[gibson_f_loc]) / GIBSON_TBE) ** 2.
        chi_squared_value += gibson_ch_sq
        # Fit to WMAP
        df += len(WMAP_F)
        wmap_ch_sq = np.sum(
            ((WMAP_TB - linear_interpolate(WMAP_F, model_frequencies, model_tb)) ** 2.) / (WMAP_TBE ** 2.)
        )
        chi_squared_value += wmap_ch_sq
        chi_squared_value = chi_squared_value[0]
    # Eliminate DoF for: saturation parameter
    df -= 1
    # Eliminate DoF for: humidity parameter
    df -= 1
    # Eliminate DoF for: offset parameter (all_others indicates sensitivity to offset; in other words NO sliding)
    df -= not all_others
    return chi_squared_value / df


def polynomial(fit, x):
    deg_proxy = len(fit) - 1
    output = np.zeros(x.size)
    for i, c in enumerate(fit):
        output += (c * (x ** float(deg_proxy - i)))
    return output


def polynomial_single(fit, x):
    deg_proxy = len(fit) - 1
    output = 0.
    for i, c in enumerate(fit):
        output += (c * (x ** float(deg_proxy - i)))
    return output


def linear_interpolate(x_targets, x_array, y_array):
    """
    Finds linear interpolation y value given a sorted list of x targets and x and y value arrays
    :param x_targets: SORTED list of x values whose y values are desired
    :param x_array: Sample x values, sorted
    :param y_array: Sample y values, corresponding to x values
    :return: y values corresponding to the list of x_targets
    """
    i = 0
    y_targets = np.zeros(len(x_targets))
    for j, x_target in enumerate(x_targets):
        while x_array[i] < x_target:
            i += 1
        if x_array[i] > x_target:
            dydx = (y_array[i] - y_array[i - 1]) / (x_array[i] - x_array[i - 1])
            dx = x_target - x_array[i - 1]
            dy = dydx * dx
            y_targets[j] = y_array[i - 1] + dy
        else:
            y_targets[j] = y_array[i]
    return y_targets


def create_parameter_space((s1, e1, step1), (s2, e2, step2), info_tuple, model_path, all_others=False):
    """
    This still works! (06/17/17, for humidity & saturation grid)
    (s1, e1, step1) are for the FIRST param
    (s2, e2 step2) are for the SECOND param
    :param info_tuple: F, TB, TBE_SLOPE, TBE_OFFSET for data
    :param model_path: String model folder name, AFTER ".../pyplanet/Output/"
    :param all_others: Whether or not this test should include WMAP, Gibson, and offset
    :return: parameter list, chi squared list
    """
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
        chi_squared_value = chi_squared(info_tuple, model_iterator.next(), all_others=all_others)
        chi_squared_values.append(chi_squared_value)
    return params, chi_squared_values


def create_parameter_space_1d(param_kws, info_tuple, model_path, all_others=False):
    """
    This has been roughly updated to work 6/22/17
    :param param_kws: looks like it supports (start, end, step) as well as a custom list
    :param info_tuple: F, TB, TBE_SLOPE, TBE_OFFSET for data
    :param model_path: String model folder name, AFTER ".../pyplanet/Output/"
    :param all_others: Whether or not this test should include WMAP, Gibson, and offset
    :return: parameter list, chi squared list
    """
    if type(param_kws) is tuple:
        params = np.arange(param_kws[0], param_kws[1], param_kws[2])
    else:
        params = param_kws
    model_path = TRUE_MODEL_PATH + model_path
    model_iterator = ModelGrid(model_path)
    chi_squared_values = []
    while model_iterator.has_next():
        chi_squared_value = chi_squared(info_tuple, model_iterator.next(), all_others=all_others)
        chi_squared_values.append(chi_squared_value)
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
        self._index = -1

    def get_index(self):
        return self._index
