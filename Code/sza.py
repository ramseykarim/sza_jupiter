"""
	==================================
	 Better version of kartos_data.py 
	==================================
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
import sys
from matplotlib.backends import backend_pdf
from scipy import constants as cst


# ----------------------------------------------------------------------------------------
def unpackcsv(csvfile, delim=','):
    myfile = open(csvfile, 'rb')
    reader = csv.reader(myfile, delimiter=delim)
    print "File: ", csvfile, " prepared for unpacking."
    return reader


"""
	Horizons
"""


# ----------------------------------------------------------------------------------------
def getdistance(row):
    l = len(row) - 2
    return float(row[l])


# ----------------------------------------------------------------------------------------
def getjdate(row):
    return float(row[1])


# ----------------------------------------------------------------------------------------
def process_horizons(reader):
    distlist = []
    datelist = []
    print "Beginning Horizons csv unpack..."
    for row in reader:
        distlist.append(getdistance(row))
        datelist.append(getjdate(row))
        print '.',
    c = 0.1202  # in AU/min
    distarr = c * np.array(distlist)
    datearr = np.array(datelist)
    print 'done'
    return datearr, distarr


"""
	SZA data from Karto
"""


# Object Oriented Attempt begins here




# ----------------------------------------------------------------------------------------
class Channel:
    def __init__(self, frequency_ghz):
        self.frequency_ghz = frequency_ghz
        self.frequency = frequency_ghz * 10. ** 9.
        self.wavelength = cst.c / self.frequency
        self.points = []

    def __str__(self):
        return str(self.frequency_ghz)

    def __repr__(self):
        return str(self)

    def methods(self):
        print "-- Methods of class Channel(frequency_ghz): --"
        print "frequency_ghz"
        print "frequency"
        print "wavelength"
        print "points"
        print "norm_distance(distances)"
        print "add_point(point)"
        print "plot_channel_vs_time(fig=11, errorbar=True, color='r')"
        print "average() = x_avg, sig_avg"
        print ""
        print ""
        print "methods()"

    def norm_distance(self, distances):
        for p, d in zip(self.points, distances):
            p.norm_distance(d)

    def add_point(self, point):
        self.points.append(point)
        point.channel = self

    def plot_channel_vs_time(self, fig=11, errorbar=True, color='r'):
        plt.figure(fig)
        point_values = getPoints()
        point_errors = getErrors()
        point_dates = np.array([p.date for p in self.points])
        if errorbar:
            plt.errorbar(point_dates, point_values, yerr=point_errors, xerr=None)
        else:
            plt.plot(point_dates, point_values, color)

    def getPoints(self):
        return np.array([p.value for p in self.points])

    def getErrors(self):
        return np.array([p.error for p in self.points])

    def average(self):
        point_values = self.getPoints()
        point_errors = self.getErrors()
        weights = 1. / point_errors
        weighted_values = point_values * weights
        x_sum = np.nansum(weighted_values)
        sig_sum = np.nansum(weights)
        x_avg = x_sum / sig_sum
        variance = np.nansum(weights * ((point_values - x_avg) ** 2.)) / sig_sum
        sig_avg = np.sqrt(variance)
        return x_avg, sig_avg


# ----------------------------------------------------------------------------------------
class Point:
    def __init__(self, value, error_percentage, date):
        if value == np.nan:
            self.value = np.nan
            self.error_percentage = np.inf
            self.error = np.inf
            self.good = False
        else:
            self.value = value
            self.error_percentage = error_percentage
            self.fix_error()
            self.good = True
        self.date = date  # MJD
        self.channel = None

    def __str__(self):
        if self.channel:
            frequencystring = " at " + str(self.channel)
        else:
            frequencystring = ""
        return "Date (MJD): " + str(self.date) + frequencystring

    def __repr__(self):
        return str(self)

    def methods(self):
        print "-- Methods of class Point(value, error_percentage, date): --"
        print "value"
        print "error_percentage"
        print "error"
        print "good"
        print "date"
        print "channel"
        print "fix_error()"
        print "norm_distance(distance)"
        print ""
        print ""
        print "methods()"

    def fix_error(self):
        self.error = self.value * self.error_percentage

    def norm_distance(self, distance):
        if self.good:
            raw_flux = self.value
            modifier = (distance / 4.04) ** 2.
            self.value = raw_flux * modifier
            self.fix_error()


# Object Oriented Attempt ends here






# ----------------------------------------------------------------------------------------
def concatstr(strlst):
    totalstr = ''
    for string in strlst:
        totalstr += str(string) + ' '
    return totalstr


# ----------------------------------------------------------------------------------------
def process_kartos(reader):
    # The text header
    header = 'Header'
    # Dates according to Karto
    dates = []
    # Frequencies
    channels = []
    # Everything
    allpoints = []
    rownum = 0
    print "Beginning SZA csv unpack..."
    for row in reader:
        if rownum == 0 or rownum == 2:
            header += '\n' + concatstr(row)
        elif rownum == 1:
            colnum = 0
            for col in row:
                if colnum % 2 == 1:
                    channels.append(float(col))
                colnum += 1
            channels.pop()
        else:
            points = []
            errors = []
            colnum = 0
            for col in row:
                if colnum == 0:
                    dates.append(float(col))
                elif colnum >= 31:
                    pass
                elif colnum % 2 == 1:
                    points.append(float(col))

                else:
                    errors.append(float(col))
                colnum += 1
            allpoints += [[points, errors]]
        rownum += 1
        print '.',
    data = np.array(allpoints)
    print 'done'
    return np.array(dates), channels, data


# ----------------------------------------------------------------------------------------
def adjust_to_object(dates, channels_list, data):
    channels = []
    print "Beginning transition from Data Cube to Object..."
    for i, ch in enumerate(channels_list):
        current_channel_object = Channel(ch)
        channel_data = data[:, :, i]
        print '.',
        for j in range(len(channel_data[:, 0])):
            value = channel_data[j, 0]
            error_percentage = channel_data[j, 1]
            date = dates[j]
            current_point_object = Point(value, error_percentage, date)
            current_channel_object.add_point(current_point_object)
        channels.append(current_channel_object)
    print "done"
    if len(channels) == 15:
        print "Success"
    else:
        print "Failed."
        print "Found only ", len(channels), " channels"
    return channels


# ----------------------------------------------------------------------------------------
def get_all_intensities(channels):
    intensities = []
    ch_errors = []
    for ch in channels:
        avg_flux, sig_avg = ch.average()
        intensities.append(avg_flux)
        ch_errors.append(sig_avg)
    intens_arr = np.array(intensities)
    error_arr = np.array(ch_errors)
    return intens_arr, error_arr


# ----------------------------------------------------------------------------------------
def normalize_all_channels(channels, distances):
    for ch in channels:
        ch.norm_distance(distances)


"""
	Plotting
"""


# def flag(times, data, distances):
#	channel_of_interest = get_channel_flux(data, 0)
#	problemindex, = np.where(channel_of_interest == np.nanmin(channel_of_interest))
#	problemindex = problemindex[0] + 1 + 4
#	print problemindex
#	return data[problemindex:, :, :], times[problemindex:], distances[problemindex:]

# ----------------------------------------------------------------------------------------
def plot_intensities(intensities, errors, chanlist, labelslist, label, xlabel="Frequency (GHz)",
                     ylabel="$T_{b}$ (K)", title="Average Intensity by Channel",
                     color='b'):
    plt.errorbar(np.array(chanlist), intensities, yerr=errors, xerr=None, fmt=".", color=color)
    #	plt.plot(chanlist, intensities, color=color)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    labelslist.append(label)


def convert_units(fluxdensities, errors, chanlist):
    flx = np.array(fluxdensities) * 10 ** (-26.)
    err = np.array(errors) * 10 ** (-26.)
    v = np.array(chanlist) * 10. ** 9.
    d_to_J_m = cst.astronomical_unit * 4.04  # m
    Req = 71492000.
    Rp = 66854000.
    SElat = 0.15 * cst.degree  # radians
    Rpp = np.sqrt(((Req) * np.sin(SElat)) ** 2. + ((Rp) * np.cos(SElat)) ** 2.)
    c = cst.c  # m/s
    h = cst.h  # Planck, SI
    pi = cst.pi
    k = cst.Boltzmann  # J/K
    Tcmb = 2.725  # Kelvins
    # Now for a bunch of math...
    first_term = 1. / (np.exp((h * v) / (k * Tcmb)) - 1)
    second_term_cst = (c ** 2. * d_to_J_m ** 2.) / (2 * h * v ** 3. * pi * Req * Rpp)
    second_term_flx = flx * second_term_cst
    second_term_err = err * second_term_cst
    inside_log_flx = 1. + 1. / (first_term + second_term_flx)
    inside_log_err = 1. + 1. / second_term_err
    Tb = ((h * v) / k) / np.log(inside_log_flx)
    Tb_err = ((h * v) / k) / np.log(inside_log_err)
    return Tb, Tb_err


# Given distance-sensitive flux densities in Jy, normalized to 4.04 AU,
# returns brightness temperature in Kelvins
# ----------------------------------------------------------------------------------------
def convert_units_old(fluxdensities, errors, chanlist, jy=True, ghz=True, neederrors=True):
    mod_flx = fluxdensities
    if neederrors:
        mod_err = errors
    # AU in meters (Wikipedia)
    AU_m = cst.astronomical_unit  # m
    d_to_J_AU = 4.04  # AU
    d_to_J_m = AU_m * d_to_J_AU  # m
    # Solid Angle calculation, given info from wikipedia
    r_J_m = 69911000.  # m
    ang_r_J = r_J_m / d_to_J_m  # rad
    solid_angle_J = np.pi * ang_r_J ** 2.  # sr
    # Convert channels to wavelength, c from Wikipedia
    c_mks = cst.c  # m s^-1
    if ghz:
        modpwr = (10. ** 9.)
    else:
        modpwr = 1.
    wavelengths = c_mks / ((np.array(chanlist)) * modpwr)  # m (= m s^-1 / (GHZ * 10^9))
    # Boltzmann constant (Wikipedia)
    kb = cst.Boltzmann  # J K^-1
    # SI conversion of flux densities and errors, both originally in Jy
    # 1 Jy = 10^-26 W m^-2 Hz^-1 = 10^-26 W s m^-2
    if jy:
        mod_flx = mod_flx / (10. ** 26.)  # W s m^-2
        if neederrors:
            mod_err = mod_err / (10. ** 26.)  # W s m^-2
    # Formula from flux density to brightness temperature:
    #   T = (wavelengths^2 * mod_flx)/(2 * kb * solid_angle_J)
    # Unit analysis:
    #   (m^2 * W s m^-2)/(J K^-1)
    #   W s K J^-1
    #   J K J^-1
    #   K
    temperature_multipliers = (wavelengths ** 2.) / (2. * kb * solid_angle_J)
    temperatures = temperature_multipliers * mod_flx
    if neederrors:
        temp_errors = temperature_multipliers * mod_err
        return temperatures, temp_errors
    else:
        return temperatures


"""
	Model Additions & Subtractions
"""


# Modeling
# ----------------------------------------------------------------------------------------
def generate_synchrotron(frequencies):  # According to J.Gibson method
    f0 = 28.5 * 10. ** 9.
    j0 = 1.5
    synchrotron = j0 * ((frequencies / f0) ** (-0.4))
    return synchrotron


def generate_cmb_old(frequencies):  # According to Planck function. Doesn't really work..
    first_term = 2. * cst.h * (frequencies ** 3.) / (cst.c ** 2.)
    second_term_inverted = np.exp((cst.h * frequencies) / (cst.Boltzmann * 2.74)) - 1.
    second_term = np.reciprocal(second_term_inverted)
    cmb = first_term * second_term
    return cmb


def generate_cmb(frequencies):
    gibson_f = [20.050, 20.295, 20.735, 20.970, 21.225, 21.550, 22.055, 22.220, 22.370, 22.945, 23.410, 23.840, 24.100,
                32.960, 40.890, 61.340, 93.820]
    gibson_cmb = [2.29, 2.28, 2.27, 2.27, 2.26, 2.26, 2.24, 2.24, 2.24, 2.23, 2.22, 2.21, 2.20, 2.0, 1.9, 1.5, 1.1]
    gibson_f = np.array(gibson_f) * 10. ** 9.
    gibson_cmb = np.array(gibson_cmb)
    deg = 2
    curve = np.polyfit(gibson_f, gibson_cmb, deg)

    def getcmb(arr):
        return curve[0] * arr ** 2 + curve[1] * arr + curve[2]

    return getcmb(frequencies)


"""
	Deal with new data (flagged & golden)
"""


# ----------------------------------------------------------------------------------------
def get_distances_times_newdata(horizons_distances, raw_times, sza_times):
    sza_distances = np.array([d for d, t in zip(horizons_distances, raw_times) if t in sza_times])
    return sza_distances


"""
	Add in Gibson KG & WMAP data!
"""


def pltkg(labelslist, label, color):
    kg_f = [20.050, 20.295, 20.735, 20.970, 21.225, 21.550, 22.055, 22.220, 22.370, 22.945, 23.410, 23.840, 24.100]
    kg_tb = [144.9, 145.6, 144.0, 143.0, 146.0, 142.2, 142.1, 137.6, 137.1, 138.6, 139.1, 141.6, 140.8]
    kg_tbe = [5.8, 5.8, 5.8, 5.7, 5.8, 5.7, 5.7, 5.5, 5.5, 5.5, 5.6, 5.7, 5.6]
    plot_intensities(kg_tb, kg_tbe, kg_f, labelslist, label, color=color)


#	plt.xscale('log', basex=10)

def pltwmap(labelslist, label):
    wmap_f = [22.44, 32.63, 40.50, 60.57, 92.71]
    wmap_tb = [134.7, 148.4, 157.1, 166.2, 174.3]
    wmap_tbe = [4.0, 2.0, 1.7, 1.5, 1.7]
    plot_intensities(wmap_tb, wmap_tbe, wmap_f, labelslist, label, color='r')


def pltwmap2(labelslist, label, color):
    wmap2_f = [22.44, 32.63, 40.50, 61.57, 92.71]
    wmap2_tb = np.array([136.5, 149.0, 157.2, 168.0, 176.0])
    wmap2_tbe = np.array([0.85, 0.77, 0.67, 0.59, 0.54])
    wmap2_tbe = wmap2_tb * wmap2_tbe / 100.
    plot_intensities(wmap2_tb, wmap2_tbe, wmap2_f, labelslist, label, color=color)


def pltgib(labelslist, label, color):
    plot_intensities([142.9], [2.3], [28.5], labelslist, label, color=color)


def pltimke(labelslist, label, color):
    imkefreq = [13.18, 14.21, 15.18, 16.21, 17.38, 8.5, 9.52, 10.46, 11.46, 4.52, 5.49, 6.5, 7.5]
    imketemp = [159.5, 154.5, 150.6, 146.4, 143.5, 187.7, 177.5, 172.4, 165.8, 247.5, 223.3, 207.6, 192.9]
    imketerr = [4.8, 4.6, 4.5, 4.4, 4.3, 5.7, 5.4, 5.2, 5., 7.4, 6.7, 6.2, 5.8]
    plot_intensities(imketemp, imketerr, imkefreq, labelslist, label, color=color)


def pltmodelold(name, color, labelslist, label):
    datafile = open(name, 'r')
    lines = datafile.readlines()
    frequencies = []
    points = []
    for line in lines:
        if not line.startswith('#'):
            pair = [float(x) for x in line.split()]
            frequencies.append(pair[0])
            points.append(pair[1])
    plt.plot(np.array(frequencies), points, color=color, linestyle='--')
    labelslist.append(label)


def pltmodel(name, color, labelslist, label):
    datafile = open(name, 'r')
    lines = datafile.readlines()
    frequencies = []
    points = []
    for line in lines:
        if not line.startswith('#'):
            pair = [float(x) for x in line.split()]
            frequencies.append(pair[0])
            points.append(pair[1])
    plt.plot(np.array(frequencies), points, color=color)
    labelslist.append(label)


# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
#                   Now doing calculations
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------






# ----------------------------------------------------------------------------------------
horiz_reader = unpackcsv('horizons_kartos/distances.csv')
time_horiz, distarr = process_horizons(horiz_reader)

jup_reader_raw = unpackcsv('jupfluxes.csv', delim=' ')
time_raw, channels_list_raw, data_raw = process_kartos(jup_reader_raw)

# jup_reader_golden = unpackcsv('jupiter_flux_golden.txt', delim=' ')
# header_golden, time_golden, channels_golden, data_golden = process_kartos(jup_reader_golden)

# jup_reader_flagged = unpackcsv('jupiter_fluxes_flagged.txt', delim=' ')
# header_flagged, time_flagged, channels_flagged, data_flagged = process_kartos(jup_reader_flagged)

jup_reader_best = unpackcsv('jupiter_fluxes_flagged2.txt', delim=' ')
time_best, channels_list, data_cube_best = process_kartos(jup_reader_best)
channels_best = adjust_to_object(time_best, channels_list, data_cube_best)

distarr_best = get_distances_times_newdata(distarr, time_raw, time_best)
normalize_all_channels(channels_best, distarr_best)

intensities_best, errors_best = get_all_intensities(channels_best)
temps_flagged2, terrs_flagged2 = convert_units(intensities_best, errors_best, channels_list)

frequencies = np.array([ch.frequency for ch in channels_best])
print "Frequencies (GHz):"
print frequencies * 10. ** (-9.)
print "Wavelengths (cm)"
wavelengths = (cst.c / frequencies) * 100.
print wavelengths
synchrotron = generate_synchrotron(frequencies)
intensities_best_s = intensities_best - synchrotron

temps_flagged2_s, terrs_flagged2_s = convert_units(intensities_best_s, errors_best, channels_list)

print "Data (K)"
print temps_flagged2_s
print "Error"
print terrs_flagged2_s

"""
fl = open('ramsey_data.txt', 'w')
fl.write("# Frequency (GHz), Wavelength (cm), T_b (K), Error (K)\n")
for f, w, t, e in zip(frequencies * 10.**(-9.), wavelengths, temps_flagged2_s, terrs_flagged2_s):
	fl.write(str(f) + ', ' + str(w) + ', ' + str(t) + ', ' + str(e) + '\n')
"""

prefix = 'models/'


def makeplots():
    labelslist = []
    pltmodel(prefix + 'pyplanetmodel_disc.dat', 'black', labelslist, "Model (Original)")
    #	pltmodel(prefix+'pyplanetmodel_80p.dat', 'red', labelslist, "Model (80% NH3 above 1 bar)")
    pltmodel(prefix + 'pyplanetmodel_50p.dat', 'purple', labelslist, "Model (50% NH3 above 1 bar)")
    #	pltmodel(prefix+'pyplanetmodel__10p.dat', 'cyan', labelslist, "Model (10% NH3 above 0.6 bar)")
    #	pltmodel(prefix+'pyplanetmodel__01p.dat', 'maroon', labelslist, "Model (1% NH3 above 0.6 bar)")
    pltmodel(prefix + 'pyplanetmodel__1pAbv06bar.dat', 'green', labelslist, "Model (1% NH3 above 0.6 bar)")
    pltmodel(prefix + 'pyplanetmodel__80pAbv2bar_10pAbv06bar.dat', 'red', labelslist,
             "Model (80% NH3 > 2bar, 10% NH3 > 0.6bar)")
    pltmodel(prefix + 'pyplanetmodel__50pAbv2bar_80pAbv06bar.dat', 'pink', labelslist,
             "Model (50% NH3 > 2bar, 80% NH3 > 0.6bar)")
    plot_intensities(temps_flagged2_s, terrs_flagged2, channels_list, labelslist, "CARMA", color='g',
                     ylabel="$T_{b}$ (K)")
    pltkg(labelslist, "KG")
    #	pltwmap(labelslist, "WMAP (Old)")
    pltwmap2(labelslist, "WMAP (Newer)")
    pltgib(labelslist, "Gibson")
    plt.legend(labelslist, loc='lower right')
    plt.xscale('log', basex=2)
    plt.xlim(19, 100)
    locs = np.arange(40) * 2 + 20
    plt.xticks(locs, [str(loc) if loc % 10 == 0 else '' for loc in locs])


#		plt.legend(["$T_{d}$", "$T_{d} + T_{CMB} - T_{synch}$"], loc='lower right')


def makeposterplots():
    labelslist = []
    #	plt.rc('font', family='serif')
    #	plt.rc('font', size=12)
    plt.rc('xtick.major', width=1)
    plt.rc('xtick.major', size=5)
    plt.rc('ytick.major', width=1)
    plt.rc('ytick.major', size=5)
    plt.rc('lines', linewidth=1.2)
    pltmodel(prefix + 'pyplanetmodel__disc.dat', 'black', labelslist, "Model (Saturated)")
    pltmodel(prefix + 'pyplanetmodel__50pAbv1bar.dat', 'purple', labelslist, "Model (50% NH3 above 1 bar)")
    pltmodel(prefix + 'pyplanetmodel__1pAbv06bar.dat', 'green', labelslist, "Model (1% NH3 above 0.6 bar)")
    pltmodel(prefix + 'pyplanetmodel__80pAbv2bar_10pAbv06bar.dat', 'blue', labelslist,
             "Model (80% NH3 > 2bar, 10% NH3 > 0.6bar)")
    pltmodel(prefix + 'pyplanetmodel__50pAbv2bar_80pAbv06bar.dat', 'maroon', labelslist,
             "Model (50% NH3 > 2bar, 80% NH3 > 0.6bar)")
    plt.errorbar(np.array(channels_list), temps_flagged2_s, yerr=terrs_flagged2, xerr=None, fmt='.', color='red')
    labelslist.append("CARMA")
    pltkg(labelslist, "KG", color='blue')
    pltwmap2(labelslist, "WMAP", color='green')
    pltgib(labelslist, "Gibson", color='purple')
    pltimke(labelslist, "de Pater", color='navy')
    plt.legend(labelslist, loc='upper right')
    plt.title("Thermal Radiation Intensity by Wavelength")
    plt.xscale('log', basex=2)
    #	print(cst.c*(1E-7)/np.array(channels_list))
    plt.xlim(4, 94)
    plt.ylim(120, 260)
    locs = np.arange(90) + 4
    plt.xticks(locs, [str(loc) if loc % 10 == 0 else '' for loc in locs])


def makeplots1():
    labelslist = []
    #	plt.rc('font', family='serif')
    #	plt.rc('font', size=12)
    plt.rc('xtick.major', width=1)
    plt.rc('xtick.major', size=5)
    plt.rc('ytick.major', width=1)
    plt.rc('ytick.major', size=5)
    plt.rc('lines', linewidth=1.2)
    pltmodel(prefix + 'pyplanetmodel__disc.dat', 'black', labelslist, "Model (Saturated)")
    pltmodel(prefix + 'pyplanetmodel__50pAbv1bar.dat', 'purple', labelslist, "Model (50% NH3 above 1 bar)")
    pltmodel(prefix + 'pyplanetmodel__1pAbv06bar.dat', 'green', labelslist, "Model (1% NH3 above 0.6 bar)")
    pltmodel(prefix + 'pyplanetmodel__80pAbv2bar_10pAbv06bar.dat', 'blue', labelslist,
             "Model (80% NH3 > 2bar, 10% NH3 > 0.6bar)")
    pltmodel(prefix + 'pyplanetmodel__50pAbv2bar_80pAbv06bar.dat', 'maroon', labelslist,
             "Model (50% NH3 > 2bar, 80% NH3 > 0.6bar)")
    plt.errorbar(np.array(channels_list), temps_flagged2_s, yerr=terrs_flagged2, xerr=None, fmt='.', color='red')
    labelslist.append("CARMA")
    pltkg(labelslist, "KG", color='blue')
    pltwmap2(labelslist, "WMAP", color='green')
    pltgib(labelslist, "Gibson", color='purple')
    pltimke(labelslist, "de Pater", color='navy')
    plt.legend(labelslist, loc='upper right')
    plt.title("Thermal Radiation Intensity by Wavelength")
    plt.xscale('log', basex=2)
    #	print(cst.c*(1E-7)/np.array(channels_list))
    plt.xlim(4, 94)
    plt.ylim(120, 260)
    locs = np.arange(90) + 4
    plt.xticks(locs, [str(loc) if loc % 10 == 0 else '' for loc in locs])


def oldmodels():
    labelslist = []
    pltmodelold(prefix + 'pyplanetmodel__disc.dat', 'black', labelslist, "Model (Saturated) - old")
    pltmodelold(prefix + 'pyplanetmodel__50pAbv1bar.dat', 'purple', labelslist, "Model (50% NH3 above 1 bar) - old")
    pltmodelold(prefix + 'pyplanetmodel__1pAbv06bar.dat', 'green', labelslist, "Model (1% NH3 above 0.6 bar) - old")
    pltmodelold(prefix + 'pyplanetmodel__80pAbv2bar_10pAbv06bar.dat', 'blue', labelslist,
                "Model (80% NH3 > 2bar, 10% NH3 > 0.6bar) - old")
    pltmodelold(prefix + 'pyplanetmodel__50pAbv2bar_80pAbv06bar.dat', 'maroon', labelslist,
                "Model (50% NH3 > 2bar, 80% NH3 > 0.6bar) - old")


def newmodels():
    labelslist = []
    #	plt.rc('font', family='serif')
    #	plt.rc('font', size=12)
    plt.rc('xtick.major', width=1)
    plt.rc('xtick.major', size=5)
    plt.rc('ytick.major', width=1)
    plt.rc('ytick.major', size=5)
    plt.rc('lines', linewidth=1.2)
    pltmodel(prefix + 'pyplanetmodel_faster_disc.dat', 'black', labelslist, "Model (Saturated)")
    pltmodel(prefix + 'pyplanetmodel_faster_50pAbv1bar.dat', 'purple', labelslist, "Model (50% NH3 above 1 bar)")
    pltmodel(prefix + 'pyplanetmodel_faster_1pAbv06bar.dat', 'green', labelslist, "Model (1% NH3 above 0.6 bar)")
    pltmodel(prefix + 'pyplanetmodel_faster_80pAbv2bar_10pAbv06bar.dat', 'blue', labelslist,
             "Model (80% NH3 > 2bar, 10% NH3 > 0.6bar)")
    pltmodel(prefix + 'pyplanetmodel_faster_50pAbv2bar_80pAbv06bar.dat', 'maroon', labelslist,
             "Model (50% NH3 > 2bar, 80% NH3 > 0.6bar)")
    plt.errorbar(np.array(channels_list), temps_flagged2_s, yerr=terrs_flagged2, xerr=None, fmt='.', color='red')
    labelslist.append("CARMA")
    pltkg(labelslist, "KG", color='blue')
    pltwmap2(labelslist, "WMAP", color='green')
    pltgib(labelslist, "Gibson", color='purple')
    pltimke(labelslist, "de Pater", color='navy')
    plt.legend(labelslist, loc='upper left')
    plt.title("Post-model upgrade")
    plt.xscale('log', basex=2)
    #	print(cst.c*(1E-7)/np.array(channels_list))
    plt.xlim(20, 95)
    plt.ylim(120, 220)
    locs = np.arange(76) + 19
    plt.xticks(locs, [str(loc) if loc % 10 == 0 else '' for loc in locs])


def modelcomp(*args):
    labelslist = []
    plt.rc('xtick.major', width=1)
    plt.rc('xtick.major', size=5)
    plt.rc('ytick.major', width=1)
    plt.rc('ytick.major', size=5)
    plt.rc('lines', linewidth=1.2)
    for i in args:
        if i == 1:
            pltmodel(prefix + 'pyplanetmodel_faster_disc.dat', 'black', labelslist, "Model (Saturated)")
            pltmodelold(prefix + 'pyplanetmodel__disc.dat', 'black', labelslist, "Model (Saturated) - old")
        elif i == 2:
            pltmodel(prefix + 'pyplanetmodel_faster_50pAbv1bar.dat', 'purple', labelslist,
                     "Model (50% NH3 above 1 bar)")
            pltmodelold(prefix + 'pyplanetmodel__50pAbv1bar.dat', 'purple', labelslist,
                        "Model (50% NH3 above 1 bar) - old")
        elif i == 3:
            pltmodel(prefix + 'pyplanetmodel_faster_1pAbv06bar.dat', 'green', labelslist,
                     "Model (1% NH3 above 0.6 bar)")
            pltmodelold(prefix + 'pyplanetmodel__1pAbv06bar.dat', 'green', labelslist,
                        "Model (1% NH3 above 0.6 bar) - old")
        elif i == 4:
            pltmodel(prefix + 'pyplanetmodel_faster_80pAbv2bar_10pAbv06bar.dat', 'blue', labelslist,
                     "Model (80% NH3 > 2bar, 10% NH3 > 0.6bar)")
            pltmodelold(prefix + 'pyplanetmodel__80pAbv2bar_10pAbv06bar.dat', 'blue', labelslist,
                        "Model (80% NH3 > 2bar, 10% NH3 > 0.6bar) - old")
        elif i == 5:
            pltmodel(prefix + 'pyplanetmodel_faster_50pAbv2bar_80pAbv06bar.dat', 'maroon', labelslist,
                     "Model (50% NH3 > 2bar, 80% NH3 > 0.6bar)")
            pltmodelold(prefix + 'pyplanetmodel__50pAbv2bar_80pAbv06bar.dat', 'maroon', labelslist,
                        "Model (50% NH3 > 2bar, 80% NH3 > 0.6bar) - old")
        elif i == 0:
            plt.errorbar(np.array(channels_list), temps_flagged2_s, yerr=terrs_flagged2, xerr=None, fmt='.',
                         color='red')
            labelslist.append("CARMA")
            pltkg(labelslist, "KG", color='blue')
            pltwmap2(labelslist, "WMAP", color='green')
            pltgib(labelslist, "Gibson", color='purple')
            pltimke(labelslist, "de Pater", color='navy')
        elif i == 6:
            pltmodel(prefix + 'pyplanetmodel__disc.dat', 'black', labelslist, "Model (Saturated) - new nh3 new h2")
            pltmodelold(prefix + 'pyplanetmodel_oldnh3newh2.dat', 'black', labelslist,
                        "Model (Saturated) - old nh3 new h2")

    plt.legend(labelslist, loc='upper left')
    plt.title("Post-model upgrade")
    plt.xscale('log', basex=2)
    plt.xlim(20, 95)
    plt.ylim(120, 220)
    locs = np.arange(76) + 19
    plt.xticks(locs, [str(loc) if loc % 10 == 0 else '' for loc in locs])


"""
plot_intensities(intensities_best, errors_best, channels_list, fig=0,
	                 ylabel="Brightness Temperature (K)", title="Average Intensity by Channel")

plt.show()
# THIS IS IMPORTANT!!!!
#time_horiz -= 2400000.5
# ^^^^^^^^^^^^^^^^^^^^^^^
# !!!!!!!!!!!!!!!!!!!!!!!!
"""

"""
plot_channels(time_raw, data_raw, channels, fig=0, title="Flagging comparison", onecolor=True, color='red')
plot_channels(time_flagged, data_flagged, channels, fig=0, onecolor=True, color='blue')
plot_channels(time_flagged2, data_flagged2, channels, fig=0, onecolor=True, color='black')
plt.plot([57109, 57109], [120, 240], 'g')
plt.xlabel("Date (MJD)")


channels_for_plotting = (2.99792458/np.array(channels))*10.
xl = "Wavelength (cm)"
plot_intensities(temps_flagged, terrs_flagged, channels_for_plotting, fig=1, ylabel="Brightness Temperature (K)",
				 color='b', xlabel=xl)
plot_intensities(temps_raw, terrs_raw, channels_for_plotting, fig=1, ylabel="Brightess Temperature (K)",
				 color='g', xlabel=xl)
plot_intensities(temps_golden, terrs_golden, channels_for_plotting, fig=1, ylabel="Brightness Temperature (K)",
				 color='r', xlabel=xl)
plot_intensities(temps_flagged2, terrs_flagged2, channels_for_plotting, fig=1, ylabel="Brightness Temperature (K)",
				 color='k', xlabel=xl)
plt.legend(['Properly Flagged', 'Flagged by Eye', '"Golden"', 'Flagged 2'], loc='upper right')


# More graphing stuff
plt.xscale('log', basex=2)
plt.ylim(130., 220.)
locs = [0.4, 0.6, 0.8, 1.0, 2.0, 4.0]
plt.xticks(locs, [str(loc) for loc in locs])
"""

"""
fl = open('freqs_improved.txt', 'w')
fr_range_line = 0.5*np.arange(30) + 20.
fr_range_further = 6.5*np.arange(10) + 35
for f in fr_range_line:
	fl.write(str(f) + '\n')
for f in fr_range_further:
	fl.write(str(f))
"""
