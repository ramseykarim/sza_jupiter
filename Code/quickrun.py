"""
Quick run script for plotting T_B data
"""
import unpack as up
import plotting as pltt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import datacontainer as dc
import models as mdl
import jackknife as jk
from itertools import cycle

MODELS = up.generate_names("../models/")


def paint_with_colors_of_wind():
    return cycle(['blue', 'green', 'red',
                  'orange', 'maroon', 'black',
                  'purple'])


def run_normal():
    u = up.Unpack().prepare().adjust()
    data = dc.DataContainer().populate()
    plt.figure(20)
    p = pltt.Plotting(wl=True)
    p.add_my_data(u)
    p.add_all_other_data()
    # p.add_model_plot("best", color='green', line_sty='--', label="45% above 1.5bar")
    p.add_model_plot("paulsolar", color='orange', line_sty='o', linewidth=2.0)
    p.add_model_plot("paulsolar_new", color='blue', line_sty='o', linewidth=2.0)
    p.add_model_plot("paulvla72x", color='cyan', line_sty='o', linewidth=2.0)
    # p.add_model_plot("psolartweaks", color='maroon', line_sty='-.')
    # p.add_model_plot("p72xtweaks", color='orange', line_sty='-.')
    comp, in_range = mdl.chi_squared(data.yield_model_comp_info(), mdl.generate_model("paulsolar"))
    comp2, in_range = mdl.chi_squared(data.yield_model_comp_info(), mdl.generate_model("paulsolar_new"))
    comp3, in_range = mdl.chi_squared(data.yield_model_comp_info(), mdl.generate_model("paulvla72x"))
    # comp4, in_range = mdl.chi_squared(data.yield_model_comp_info(), mdl.generate_model("50abv2bar_new"))
    # u.print_points()
    # print comp
    # print comp2, comp3, comp4
    # plt.xlim([20, 100])
    plt.figure(20)
    plt.ylim([120, 200])
    p.show_prep()
    p.show()


def run_jackknife():
    data = dc.DataContainer().populate()
    f = np.concatenate([data.frequencies, data.frequencies])
    df = 0.05
    f -= df * len(jk.JACKKNIFE_TESTS)
    u = jk.JackknifeUnpack()
    plt.figure(10)
    plt.title("Jackknife Tests")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Variation from $T_{b}$ (K)")
    legend = []
    max_error = np.max(data.tbe_offset)
    plt.axhline(0, linestyle='--', color="black", label="_nolegend_")
    # plt.plot(data.frequencies, data.tb, color="black", label='_nolegend_')
    for t in jk.JACKKNIFE_TESTS:
        u.accept_test(jk.JACKKNIFE_TESTS[t][0]).adjust()
        tb, unused, tbe_offset = u.get_temperatures()
        max_error = max(max_error, np.max(tbe_offset))
        tb -= data.tb
        tb = np.concatenate([tb, tb])
        tbe_offset = np.concatenate([data.tbe_offset, tbe_offset])
        plt.figure(10)
        plt.errorbar(f, tb, tbe_offset,
                     fmt=jk.JACKKNIFE_TESTS[t][1],
                     color=jk.JACKKNIFE_TESTS[t][2])
        legend.append(t)
        f += df
    plt.figure(10)
    b = jk.BootstrapAverageContainer(f[:f.size / 2], u)
    b.cycle(data.tb)
    plt.plot(b.f, b.t, ',', color='k')
    legend = ["Bootstrap"] + legend
    # plt.ylim([0 - 2 * max_error, 0 + 2 * max_error])
    plt.legend(legend)
    plt.show()


def run_chisq():
    data = dc.DataContainer().populate()
    model_dict = {}
    for m in MODELS:
        # comp = mdl.old_chi_sq(data.old_yield_model_comp(m))
        comp, in_range = mdl.chi_squared(data.yield_model_comp_info(), mdl.generate_model(m))
        model_dict[m] = comp
    for (key, value) in [(key, value) for (key, value) in sorted(model_dict.items(), key=lambda x: x[1])]:
        if value:
            print key + " - > ", value


def run_models():
    p = pltt.Plotting()
    plt.figure()
    u = up.Unpack().prepare().adjust()
    for m in MODELS:
        p.add_model_plot(m, color="default")
    plt.title("models")
    p.add_my_data(u)
    p.add_all_other_data()
    p.show_prep(location='lower left')
    p.show()


def run_raw():
    """
    Figure-ready, raw data before and after normalization
    """
    u = up.Unpack().prepare()
    plt.figure()
    plt.subplot(121)
    plt.title("Raw Data - Before Normalization")
    plt.xlabel("Date (MJD)")
    plt.ylabel("Flux density (jy)")
    date_max = np.max(u.dates_copss_array)
    date_min = np.min(u.dates_copss_array)
    date_range = date_max - date_min
    l_padding = 0.05
    r_padding = 0.2
    for ch in u.channel_obj_list:
        plt.errorbar(u.dates_copss_array, ch.flux, yerr=ch.error, fmt='.')
    plt.xlim(date_min - l_padding * date_range, np.max(u.dates_copss_array) + l_padding * date_range)
    plt.subplot(122)
    plt.title("Raw Data - After Normalization")
    plt.xlabel("Date (MJD)")
    plt.ylabel("Flux density (jy)")
    u.adjust()
    for ch in u.channel_obj_list:
        plt.errorbar(u.dates_copss_array, ch.flux, yerr=ch.error, fmt='.')
        plt.text(date_max + 3,
                 np.nanmean(ch.flux),
                 "%.3f GHz" % ch.frequency_ghz, fontsize=8)
    plt.xlim(date_min - l_padding * date_range, np.max(u.dates_copss_array) + r_padding * date_range)
    plt.show()


def run_model_grid_old_gas_profile():
    """
    This is really only good to get an idea of what the models look like in Pmult/Pcutoff space
    :return nothing
    """
    markers = cycle(['D', 'o', '^'])
    param1 = [float(x) / 100. for x in range(70, 100, 10)]
    param2 = [float(x) / 2. for x in range(1, 7, 1)]
    params = []
    # u = up.Unpack().prepare().adjust()
    # p = pltt.Plotting()
    # p.add_my_data(u)
    # p.add_all_other_data()
    # p.add_model_plot("saturated", color='orange', line_sty='-.')
    #
    # debug_model_path = "/home/ramsey/Documents/Research/pyplanet/Output/Jupiter_Spectrum_20161130_1413.dat"
    # debug_data = np.loadtxt(debug_model_path)
    # f, x = debug_data[:, 0], debug_data[:, 1]
    # plt.plot(f, x, '--', color='maroon', label="_nolabel_")
    #
    # p.show_prep()
    for p1 in param1:
        for p2 in param2:
            params.append((p1, p2))
    model_path = mdl.TRUE_MODEL_PATH + "Jupiter_frequency_20161130_0533/Jupiter000"
    suffix = ".dat"
    names = ["%02d" % x for x in range(18)]
    legend = []
    m = None
    df = 0.2 / len(names)
    for i, name in enumerate(names):
        if i % 7 == 0:
            m = markers.next()
        full_name = model_path + name + suffix
        data = np.loadtxt(full_name)
        plt.plot(data[:, 0] - 0.1 + (df * i), data[:, 1], m)
        legend.append("(P*=%.1f, P<%.1f)" % (params[i][0], params[i][1]))
    plt.legend(legend)
    plt.show()


def run_model_grid():
    #  Generate parameter space as it was during this run (12/01/2016) # 1
    data = dc.DataContainer().populate()
    data_model_info = data.yield_model_comp_info()
    params_1, chi_squared_values_1 = mdl.create_parameter_space((0.6, 0.9, 0.05), (0.5, 1.5, 0.1),
                                                                data_model_info, "Jupiter_frequency_20161201_1905/")
    params_2, chi_squared_values_2 = mdl.create_parameter_space((0.3, 0.6, 0.05), (0.5, 1.5, 0.1),
                                                                data_model_info, "Jupiter_frequency_20161201_2334/")
    params_3, chi_squared_values_3 = mdl.create_parameter_space((0.9, 1.3, 0.05), (0.5, 1.5, 0.1),
                                                                data_model_info, "Jupiter_frequency_20161204_1635/")
    params_4, chi_squared_values_4 = mdl.create_parameter_space((0.3, 0.6, 0.05), (1.5, 2.5, 0.1),
                                                                data_model_info, "Jupiter_frequency_20161205_1007/")
    params_5, chi_squared_values_5 = mdl.create_parameter_space((0.3, 1.3, 0.05), (0.0, 0.5, 0.1),
                                                                data_model_info, "Jupiter_frequency_20161213_1403/")
    params_6, chi_squared_values_6 = mdl.create_parameter_space((1.3, 1.5, 0.05), (0.5, 1.0, 0.1),
                                                                data_model_info, "Jupiter_frequency_20161213_1616/")
    params = params_1 + params_2 + params_3 + params_4 + params_5 + params_6
    chi_squared_values = chi_squared_values_1 + chi_squared_values_2 + \
                         chi_squared_values_3 + chi_squared_values_4 + \
                         chi_squared_values_5 + chi_squared_values_6

    # Plot up results
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    saturation = [x[0] for x in params]
    cutoff = [x[1] for x in params]
    chi_squared_value = [x[0] for x in chi_squared_values]
    test_array = np.array(chi_squared_value)
    lowest_val = np.where(chi_squared_value == np.min(chi_squared_value))
    test_array1 = [params[i] for i in lowest_val[0]]
    print "Lowest Chi Squared Values:"
    for (p1, p2) in test_array1:
        print "Sat: %.2f, Cutoff:%.2f" % (round(p1, 2), round(p2, 1))
    in_range = [x[1] for x in chi_squared_values]
    ax.plot([x for x, i in zip(saturation, in_range) if i],
            [x for x, i in zip(cutoff, in_range) if i],
            zs=[x for x, i in zip(chi_squared_value, in_range) if i],
            color='black',
            marker='o', linestyle='None')
    ax.plot([x for x, i in zip(saturation, in_range) if not i],
            [x for x, i in zip(cutoff, in_range) if not i],
            zs=[x for x, i in zip(chi_squared_value, in_range) if not i],
            color='red',
            marker='o', linestyle='None')

    plt.xlabel("Saturation (multiplier)")
    plt.ylabel("Cutoff (bar)")

    #  Visually compare to points
    # plt.figure()
    # model_iterator_1.reinitialize()
    # while model_iterator_1.has_next():
    #     f, tb = model_iterator_1.next()
    #     plt.plot(f, tb)
    # p = pltt.Plotting()
    # p.add_my_data(data)
    plt.show()


def run_air_mass_investigation():
    u = up.Unpack().prepare().adjust()
    plt.figure(10)
    colors = paint_with_colors_of_wind()

    #  Flux vs elevation
    j_el_raw = u.get_j_el()
    freqs = []
    slopes = []
    for ch in u.channel_obj_list:
        freqs.append(ch.frequency_ghz)
        flux = ch.flux[np.where(ch.flux < np.inf)]
        j_el = j_el_raw[np.where(ch.flux < np.inf)]
        err = ch.error[np.where(ch.flux < np.inf)]
        color = colors.next()
        plt.errorbar(j_el, flux, yerr=err, fmt='.', color=color)
        e_fit = np.polyfit(j_el, flux, deg=1)
        e_fit_line = j_el * e_fit[0] + (j_el ** 0) * e_fit[1]
        slopes.append(e_fit[0])
        plt.plot(j_el, e_fit_line, '--', color=color)
        plt.text(72, np.mean(flux), "%.2f GHz" % ch.frequency_ghz, fontsize=8)
    plt.title("Flux vs Jupiter Elevation (obs date & JPL Horizons) and linear fits")
    plt.xlabel("Elevation (degrees)")
    plt.ylabel("Flux (jy) (distance adjusted)")

    #  Elevation vs time
    plt.figure(11)
    dates = np.array(u.correlated_dates)
    plt.plot(dates, j_el_raw, '.', color="red")
    fit = np.polyfit(dates, j_el_raw, deg=1)
    fit_line = dates * fit[0] + (dates ** 0. * fit[1])
    plt.plot(dates, fit_line, '--', color="black")
    plt.ylim([0, 90])
    plt.title("Elevation vs Time and linear fit")
    plt.xlabel("Time (JD)")
    plt.ylabel("Elevation (degrees)")

    #  Flux/elevation fit slope vs frequency
    plt.figure(12)
    freqs = np.array(freqs)
    slopes = np.array(slopes)
    s_fit = np.polyfit(freqs, slopes, deg=1)
    s_fit_line = freqs * s_fit[0] + (freqs ** 0.) * s_fit[1]
    plt.plot(freqs, slopes, 'x', color='maroon')
    plt.plot(freqs, s_fit_line, '--', color='cyan')
    plt.ylim([-.1, .2])
    # noinspection PyStringFormat
    plt.text(30, 0, "Linear slope: %.4f jy/degree/GHz" % s_fit[0])
    plt.title("Flux/Elevation Linear Fit Slope vs Frequency and linear fit")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Flux-Elevation fit slope (jy/degree)")

    plt.show()


run_model_grid()
