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
from scipy.stats.mstats import gmean

MODELS = up.generate_names("../models/")


def paint_with_colors_of_wind():
    return cycle(['blue', 'green', 'red',
                  'orange', 'maroon', 'black',
                  'purple'])


def run_normal():
    u = up.Unpack().prepare().adjust()
    wl = True
    p = pltt.Plotting(wl=wl)
    # Playing with P<0.5bar
    # p.add_model_plot("0d5_40pct_1d6", color='gold', line_sty='--', label="0.5<(40%x)<1.6 bar")
    # p.add_model_plot("0d5_40pct_1d6_60pct", color='cyan', line_sty='--', label="0.5<(40%x)<1.6<(60%x)")
    # p.add_model_plot("0d5_30pct_1d6_40pct", color='blue', line_sty='--', label="0.5<(30%x)<1.6<(40%x)")
    # p.add_model_plot("0d5_60pct_2d6_40pct", color='green', line_sty='--', label="0.5<(60%x)<2.6<(40%x)")
    p.add_model_plot("0d5_60pct", color='cyan', line_sty='-', label="0.5<(60%x)")
    p.add_model_plot("0d5_50pct", color='blue', line_sty='-', label="0.5<(50%x)")
    p.add_model_plot("0d5_45pct", color='green', line_sty='-', label="0.5<(45%x)")
    p.add_model_plot("0d5_40pct", color='orange', line_sty='-', label="0.5<(40%x)")
    p.add_model_plot("0d5_30pct", color='maroon', line_sty='-', label="0.5<(30%x)")
    # p.add_model_plot("0d5_45pct_1d6_50pct", color='magenta', line_sty='--', label="0.5<(45%x)<1.6<(50%x)")
    # p.add_model_plot("0d5_50pct_0d8_45pct", color='red', line_sty='--', label="0.5<(50%x)<0.8<(45%x)")
    # p.add_model_plot("0d5_50pct_1d0_45pct", color='blue', line_sty='--', label="0.5<(50%x)<1.0<(45%x)")
    # p.add_model_plot("0d5_50pct_1d2_45pct", color='green', line_sty='--', label="0.5<(50%x)<1.2<(45%x)")
    # p.add_model_plot("0d5_50pct_1d4_45pct", color='black', line_sty='--', label="0.5<(50%x)<1.4<(45%x)")
    # p.add_model_plot("0d5_50pct_1d6_45pct", color='red', line_sty='--', label="0.5<(50%x)<1.6<(45%x)")
    # p.add_model_plot("0d5_60pct_1d6_40pct", color='red', line_sty='--', label="0.5<(60%x)<1.6<(40%x)")

    # p.add_model_plot("10pct_0d5_40pct_1d6", color='cyan', line_sty='--', label="\" & 10%x Nominal P<0.5 bar")
    # p.add_model_plot("200pct_0d5_40pct_1d6", color='red', line_sty='--', label="\" & 200%x Nominal P<0.5 bar")

    # Who knows anymore
    # p.add_model_plot("paulsolar_old", color='blue', line_sty='o', linewidth=2.0)
    # p.add_model_plot("paulsolar", color='orange', line_sty='-', linewidth=2.0, label="1x Solar abundance")
    # p.add_model_plot("paulvla72x", color='blue', line_sty='-', linewidth=2.0, label="Nominal abundance")
    # p.add_model_plot("0d5_45pct_1d8", color='blue', line_sty='-', linewidth=3.0, label="45%x Nominal 0.5<P<1.8 bar")

    p.add_my_data(u)
    p.add_all_other_data()

    plt.ylim([120, 200])

    if wl:
        plt.title("Brightness Temperature (disk averaged) by Wavelength")
        plt.xlabel("$\lambda$ (cm)")
    else:
        plt.title("Brightness Temperature (disk averaged) by Frequency")
        plt.xlabel("$f$ (GHz)")
        plt.xlim([20, 100])

    plt.ylabel("$T_{b}$ (K)")
    p.show_prep()
    p.show()


def run_solar_vla():
    p = pltt.Plotting(wl=True)
    p.add_model_plot("paulsolar", color='orange', line_sty='-', linewidth=2.0)
    p.add_model_plot("paulvla72x", color='blue', line_sty='-', linewidth=2.0)
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
        f_arr = [np.mean(tb)]
        tb = [np.mean(tb)]
        tbe_offset = [np.mean(tbe_offset)]
        plt.errorbar(f_arr, tb, tbe_offset,  # FIXME poorly done please fix!!
                     fmt=jk.JACKKNIFE_TESTS[t][1],
                     color=jk.JACKKNIFE_TESTS[t][2])
        legend.append(t)
        f += df
    plt.figure(10)
    # b = jk.BootstrapAverageContainer(f[:f.size / 2], u)
    # b.cycle(data.tb)
    # plt.plot(b.f, b.t, ',', color='k')
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

    # w/o WMAP
    # params_1, chi_squared_values_1 = mdl.create_parameter_space((0.6, 0.9, 0.05), (0.5, 1.5, 0.1),
    #                                                             data_model_info, "mgrid_noWMAP_1/")
    # params_2, chi_squared_values_2 = mdl.create_parameter_space((0.3, 0.6, 0.05), (0.5, 1.5, 0.1),
    #                                                             data_model_info, "mgrid_noWMAP_2/")
    # params_3, chi_squared_values_3 = mdl.create_parameter_space((0.9, 1.3, 0.05), (0.5, 1.5, 0.1),
    #                                                             data_model_info, "mgrid_noWMAP_3/")
    # params_4, chi_squared_values_4 = mdl.create_parameter_space((0.3, 0.6, 0.05), (1.5, 2.5, 0.1),
    #                                                             data_model_info, "mgrid_noWMAP_4/")
    # params_5, chi_squared_values_5 = mdl.create_parameter_space((0.3, 1.3, 0.05), (0.0, 0.5, 0.1),
    #                                                             data_model_info, "mgrid_noWMAP_5/")
    # params_6, chi_squared_values_6 = mdl.create_parameter_space((1.3, 1.5, 0.05), (0.5, 1.0, 0.1),
    #                                                             data_model_info, "mgrid_noWMAP_6/")
    # params = params_1 + params_2 + params_3 + params_4 + params_5 + params_6
    # chi_squared_values = chi_squared_values_1 + chi_squared_values_2 + chi_squared_values_3 + \
    #                      chi_squared_values_4 + chi_squared_values_5 + chi_squared_values_6

    # With WMAP
    w = True
    params_1, chi_squared_values_1 = mdl.create_parameter_space((0.3, 1.1, 0.05), (0.6, 1.2, 0.1),
                                                                data_model_info, "mgrid_WMAP_1/",
                                                                w=w)
    params_2, chi_squared_values_2 = mdl.create_parameter_space((0.3, 0.6, 0.05), (1.2, 1.8, 0.1),
                                                                data_model_info, "mgrid_WMAP_2/",
                                                                w=w)
    params_3, chi_squared_values_3 = mdl.create_parameter_space((0.35, 0.5, 0.05), (1.8, 2.1, 0.1),
                                                                data_model_info, "mgrid_WMAP_3/",
                                                                w=w)

    params = params_1 + params_2 + params_3
    chi_squared_values = chi_squared_values_1 + chi_squared_values_2 + chi_squared_values_3

    # Plot up results
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    saturation = [x[0] for x in params]
    cutoff = [x[1] for x in params]
    chi_squared_value = [x[0] for x in chi_squared_values]
    test_array = np.array(chi_squared_value)
    lowest_val = np.where(test_array == np.min(test_array))
    test_array1 = [params[i] for i in lowest_val[0]]
    print "Chi Sq Min = %.1f" % round(np.min(test_array))
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

    plt.figure()
    plt.plot(test_array)

    #  Visually compare to points
    # plt.figure()
    # model_iterator_1.reinitialize()
    # while model_iterator_1.has_next():
    #     f, tb = model_iterator_1.next()
    #     plt.plot(f, tb)
    # p = pltt.Plotting()
    # p.add_my_data(data)
    plt.show()


def run_model_grid_1d():
    """
    Abandons P_cutoff and focuses on the 1D case of abundance multiplier only.
    :return:
    """
    data = dc.DataContainer().populate()
    data_model_info = data.yield_model_comp_info()
    deg = 7  # To get a good interpolation of model fit
    color_wheel = paint_with_colors_of_wind()
    markers = ['o', 's', 'v', '^', '*']

    plt.figure()
    fits = []
    fit_errors = [[], []]
    w = True
    offset = True
    labels = ['relative only', 'absolute + WMAP + Gibson']
    x_lim = [0, 1]
    for i in range(2):
        w = not w
        offset = not offset
        # offset = not offset if not w else offset  # Switches every 2 turns (2x2 boolean grid with w)
        # Either These ###########
        params_1, chi_squared_values_1 = mdl.create_parameter_space_1d((0.25, 0.75, 0.05),
                                                                   data_model_info, "mgrid_abundance_only_1/",
                                                                   w=w, offset=offset)  # 0.5 < P < 1.8
        param0 = np.array([0.001, 0.01, 0.03], dtype=float)
        param2 = np.arange(0.05, 0.25, 0.05)  # Saturation coefficient
        param1 = np.arange(0.25, 0.75, 0.05)  # not used
        param3 = np.arange(0.75, 1.15, 0.05)
        params_final = np.concatenate([param0, param2, param3])
        params_2, chi_squared_values_2 = mdl.create_parameter_space_1d(params_final,
                                                                       data_model_info, "mgrid_abundance_only_2/",
                                                                       w=w, offset=offset)  # 0.5 < P < 1.8
        params_all = np.concatenate([params_1, params_2])
        chi_squared_value_all = np.concatenate([chi_squared_values_1, chi_squared_values_2])

        # Or This ############
        # param0 = np.array([0.001, 0.01, 0.03], dtype=float)
        # param2 = np.arange(0.05, 0.25, 0.05)  # Saturation coefficient
        # params_1 = np.arange(0.25, 0.75, 0.05)  # not used
        # param3 = np.arange(0.75, 1.15, 0.05)
        # params_final = np.concatenate([param0, param2, params_1, param3])
        # params_all, chi_squared_value_all = mdl.create_parameter_space_1d(params_final,
        #                                                                   data_model_info, "mgrid_abundance_nocutoff_1/",
        #                                                                   w=w, offset=offset)  # 0.5 < P < 1.8
        # chi_squared_values_1 = chi_squared_value_all[np.arange(params_all.shape[0])[np.in1d(params_all, params_1)]]

        clr = color_wheel.next()
        saturation = np.array(params_1)
        chi_squared_value = chi_squared_values_1
        # plt.subplot(211)
        plt.plot(params_all, chi_squared_value_all, markers[i], color=clr)
        fit = np.polyfit(saturation, chi_squared_value, deg=deg)
        x_lim = [0.2, 0.8]
        x_range = np.arange(x_lim[0], x_lim[1], 0.01)
        output = np.zeros(len(x_range))
        for j, coefficient in enumerate(fit):
            output += coefficient * (x_range ** (deg - j))
        plt.plot(x_range, output, '--', color=clr)
        where_min = np.argmin(output)
        fits.append(x_range[where_min])
        # noinspection PyTypeChecker
        where_lo_error = np.where(
            np.min((output[where_min] * 2. - output[:where_min]) ** 2) == (
                output[where_min] * 2. - output[:where_min]) ** 2)
        fit_errors[0].append(fits[i] - x_range[:where_min][where_lo_error])
        # noinspection PyTypeChecker
        where_hi_error = np.where(
            np.min((output[where_min:] - output[where_min] * 2.) ** 2) == (
                output[where_min:] - output[where_min] * 2.) ** 2)
        fit_errors[1].append(x_range[where_min:][where_hi_error] - fits[i])
        plt.errorbar([fits[i]], [output[where_min] * 3.], xerr=[fit_errors[0][i], fit_errors[1][i]], fmt=markers[i], color=clr, label=labels[i])
        # plt.subplot(212)
        # plt.errorbar([fits[i]], [i], xerr=[fit_errors[0][i], fit_errors[1][i]], fmt=markers[i], color=clr, label=label)
        # print fits[i]

    print '-'
    # plt.subplot(211)
    x_lim = [0, 1]
    plt.xlim(x_lim)
    plt.xlabel("Saturation Coefficient")
    plt.ylabel("Reduced $\chi^{2}$")
    plt.title("Model comparison to CARMA and WMAP data")
    plt.yscale('log')
    # plt.subplot(212)
    plt.legend(loc='lower right')
    # plt.ylim([-1, 4])
    # plt.xlim(x_lim)
    plt.show()


def run_air_mass_investigation():
    u = up.Unpack().prepare().adjust()
    colors = paint_with_colors_of_wind()

    #  Flux vs elevation
    j_el_raw = u.get_j_el()
    freqs = []
    slopes = []
    flat_slopes = []
    for ch in u.channel_obj_list:
        plt.figure(10)
        plt.subplot(121)
        freqs.append(ch.frequency_ghz)
        flux = ch.flux[np.where(ch.flux < np.inf)]
        j_el = j_el_raw[np.where(ch.flux < np.inf)]
        err = ch.error[np.where(ch.flux < np.inf)]
        color = colors.next()
        plt.errorbar(j_el, flux, yerr=err, fmt='.', color=color)
        t0 = flux[flux.size - 1]
        e0 = j_el[j_el.size - 1]
        e_fit = np.polyfit(j_el, flux, deg=1)
        e_fit_line = j_el * e_fit[0] + (j_el ** 0) * e_fit[1]
        factors = e_fit[0] * (j_el - e0) / t0 + 1
        flux_adj = flux / factors
        slopes.append(e_fit[0])
        plt.plot(j_el, e_fit_line, '-', color=color)
        plt.text(72, np.mean(flux), "%.2f GHz" % ch.frequency_ghz, fontsize=8)
        # plt.figure(13)
        # plt.errorbar(j_el, flux_adj, yerr=err, fmt='.', color=color)
        # plt.text(72, np.mean(flux), "%.2f GHz" % ch.frequency_ghz, fontsize=8)
        # e_fit_flat = np.polyfit(j_el, flux_adj, deg=1)
        # e_fit_flat_line = e_fit_flat[0] * j_el + (j_el ** 0.) * e_fit_flat[1]
        # flat_slopes.append(e_fit_flat[0])
        # plt.plot(j_el, e_fit_flat_line, '--', color=color)
    plt.figure(10)
    plt.title("Flux vs Jupiter Elevation and Linear Fits")
    plt.xlabel("Elevation (degrees)")
    plt.ylabel("Normalized Flux (jy)")
    # plt.figure(13)
    # plt.title("Flux vs Jupiter Elevation (obs date & JPL Horizons) flattened")
    # plt.xlabel("Elevation (degrees)")
    # plt.ylabel("Flux (jy) (distance adjusted)")

    #  Elevation vs time
    # plt.figure(11)
    # dates = np.array(u.correlated_dates)
    # plt.plot(dates, j_el_raw, '.', color="red")
    # fit = np.polyfit(dates, j_el_raw, deg=1)
    # fit_line = dates * fit[0] + (dates ** 0. * fit[1])
    # plt.plot(dates, fit_line, '--', color="black")
    # plt.ylim([0, 90])
    # plt.title("Elevation vs Time and linear fit")
    # plt.xlabel("Time (JD)")
    # plt.ylabel("Elevation (degrees)")

    #  Flux/elevation fit slope vs frequency
    plt.figure(10)
    plt.subplot(122)
    freqs = np.array(freqs)
    slopes = np.array(slopes)
    s_fit = np.polyfit(freqs, slopes, deg=1)
    s_fit_line = freqs * s_fit[0] + (freqs ** 0.) * s_fit[1]
    plt.plot(freqs, slopes, 'D', color='maroon', markersize=5.0)
    plt.plot(freqs, s_fit_line, '--', color='blue', linewidth=2.0)
    plt.ylim([-.1, .2])
    # noinspection PyStringFormat
    plt.text(30, 0, "Linear slope: %.4f jy/degree/GHz" % s_fit[0])
    plt.title("Flux/Elevation Linear Fit Slope vs Frequency and Linear Fit")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Flux-Elevation Fit Slope (jy/degree)")

    # plt.figure(14)
    # flat_slopes = np.array(flat_slopes)
    # s_flat_fit = np.polyfit(freqs, flat_slopes, deg=1)
    # s_flat_fit_line = freqs * s_flat_fit[0] + (freqs ** 0.) * s_flat_fit[1]
    # plt.plot(freqs, flat_slopes, 'x', color='maroon')
    # plt.plot(freqs, s_flat_fit_line, '--', color='cyan')
    # plt.ylim([-.1, .2])
    # # noinspection PyStringFormat
    # plt.text(30, 0, "Linear slope: %.4f jy/degree/GHz" % s_flat_fit[0])
    # plt.title("Flux/Elevation Linear Fit Slope vs Frequency and linear fit, fixed")
    # plt.xlabel("Frequency (GHz)")
    # plt.ylabel("Flux-Elevation fit slope (jy/degree)")

    plt.show()


def run_atmosphere_graphs():
    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulvla72x',
                         skip_header=1)
    print data.shape
    plt.figure()
    layer = data[:, 0]
    t = data[:, 1]
    p = data[:, 2]
    h2 = data[:, 3]
    he = data[:, 4]
    ch4 = data[:, 5]
    nh3 = data[:, 6]
    h2o = data[:, 7]
    plt.subplot(121)
    plt.title("Temperature-Pressure Profile")

    # TP Graph
    plt.plot(t, p)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1.e2, 1.e-2])
    plt.xlim(([1.e2, 1.e3]))
    plt.xlabel("T (K)")
    plt.ylabel("P (bar)")

    plt.subplot(122)
    plt.title("NH3 Abundance Models")

    #  NOMINAL ABUNDANCE
    # plt.plot(h2, p)
    # plt.plot(he, p)
    # plt.plot(ch4, p)
    plt.plot(nh3, p, label='NH3 $-$ Nominal')
    print nh3[p - 1 == np.min(np.abs(p - 1))]
    # plt.plot(h2o, p)
    # plt.legend(['h2', 'he', 'ch4', 'nh3', 'h2o'])
    tru_nh3 = np.array(nh3)

    # MY FUCKED UP THREE STAGE ABUNDANCE
    # plim = np.where(p < 0.5)
    # # nh3[plim] *= 2.0
    # # plt.plot(nh3[plim], p[plim], '--', label='__nolabel__')
    # plim = np.where((p >= 0.5) & (p < 0.8))
    # nh3[plim] *= 0.5
    # # plt.plot(nh3[plim], p[plim], '--', label='__nolabel__')
    # # for u, v in zip(p[plim], nh3[plim]):
    # #     print u, 'bar: ', v
    # plim = np.where(p >= 0.8)
    # nh3[plim] *= 0.45
    # # plt.plot(nh3[plim], p[plim], '--', label='__nolabel__')
    # plim = np.where((p >= 0.8) & (p < 2.0))
    # # for u, v in zip(p[plim], nh3[plim]):
    # #     print u, 'bar: ', v
    # plt.plot(nh3, p, '--', label='50%/45% (.5/.8bar)')

    print '=='
    #  60 PCT
    nh3 = np.copy(tru_nh3)
    print nh3[-1]
    print nh3[0]
    plim = np.where(p >= 0.5)
    nh3[plim] *= 0.6
    print '-'
    print nh3[0]
    print nh3[p - 1 == np.min(np.abs(p - 1))]
    print nh3[-1]
    print '--'
    plt.plot(nh3, p, '--', label='_nolegend_')
    plt.text(nh3[p == np.max(p)], 10*0.6, "60%")

    #  50 PCT
    nh3 = np.copy(tru_nh3)
    print nh3[-1]
    print nh3[0]
    plim = np.where(p >= 0.5)
    nh3[plim] *= 0.5
    print '-'
    print nh3[0]
    print nh3[p - 1 == np.min(np.abs(p - 1))]
    print nh3[-1]
    print '--'
    plt.plot(nh3, p, '--', label='_nolegend_')
    plt.text(nh3[p == np.max(p)], 10*0.5, "50%")

    #  45 PCT
    nh3 = np.copy(tru_nh3)
    print nh3[-1]
    print nh3[0]
    plim = np.where(p >= 0.5)
    nh3[plim] *= 0.45
    print '-'
    print nh3[0]
    print nh3[p - 1 == np.min(np.abs(p - 1))]
    print nh3[-1]
    print '--'
    plt.plot(nh3, p, '--', label='_nolegend_')
    plt.text(nh3[p == np.max(p)], 10*0.45, "45%")

    #  40 PCT
    nh3 = np.copy(tru_nh3)
    print nh3[-1]
    print nh3[0]
    plim = np.where(p >= 0.5)
    nh3[plim] *= 0.4
    print '-'
    print nh3[0]
    print nh3[p - 1 == np.min(np.abs(p - 1))]
    print nh3[-1]
    print '--'
    plt.plot(nh3, p, '--', label='_nolegend_')
    plt.text(nh3[p == np.max(p)], 10*0.4, "40%")

    #  30 PCT
    nh3 = np.copy(tru_nh3)
    print nh3[-1]
    print nh3[0]
    plim = np.where(p >= 0.5)
    nh3[plim] *= 0.3
    print '-'
    print nh3[0]
    print nh3[p - 1 == np.min(np.abs(p - 1))]
    print nh3[-1]
    print '--'
    plt.plot(nh3, p, '--', label='_nolegend_')
    plt.text(nh3[p == np.max(p)], 10*0.3, "30%")

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.e-8, 1.e-2])
    plt.ylim([1.e2, 1.e-2])
    plt.xlabel("Fractional Abundance")

    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulSolar_asplund_ws',
                         skip_header=1)
    p = data[:, 2]
    nh3 = data[:, 6]
    plt.plot(nh3, p, label='Solar')

    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulclvla72x',
                         skip_header=1)
    # plt.subplot(122)

    layer = data[:, 0]
    t = data[:, 1]
    p = data[:, 2]
    soln = data[:, 3]
    h2o = data[:, 4]
    nh4sh = data[:, 5]
    nh3 = data[:, 6]

    def get_x_pos(x, y):
        # return x[x == np.max(x)]
        return gmean([x[x == np.max(x)], np.array([plt.xlim()[0]])])

    def get_y_pos(x, y):
        return np.sum(y[x != 0] * x[x != 0])/np.sum(x[x != 0])

    plt.plot(nh3, p, label='_nolegend_')
    plt.text(get_x_pos(nh3, p), get_y_pos(nh3, p), "NH3")
    plt.plot(soln, p, label='_nolegend_')
    plt.text(get_x_pos(soln, p), get_y_pos(soln, p), "SOLN")
    plt.plot(h2o, p, label='_nolegend_')
    plt.text(get_x_pos(h2o, p), get_y_pos(h2o, p), "H2O")
    plt.plot(nh4sh, p, label='_nolegend_')
    plt.text(get_x_pos(nh4sh, p), get_y_pos(nh4sh, p), "NH4SH")
    # plt.xlabel("Cloud (?)")
    # plt.ylabel("Pressure (bar)")
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.ylim([1.e2, 1.e-2])
    plt.legend()

    plt.show()


def run_atmosphere_graphs_2():
    """
    Trying to create a Jupiter profile from proto-solar value
    :return: None
    """
    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulvla72x',
                         skip_header=1)
    print data.shape
    plt.figure()
    layer = data[:, 0]
    t = data[:, 1]
    p = data[:, 2]
    h2 = data[:, 3]
    he = data[:, 4]
    ch4 = data[:, 5]
    nh3 = data[:, 6]
    h2o = data[:, 7]
    plt.subplot(121)
    plt.title("Temperature-Pressure Profile")

    # TP Graph
    plt.plot(t, p)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1.e2, 1.e-2])
    plt.xlim(([1.e2, 1.e3]))
    plt.xlabel("T (K)")
    plt.ylabel("P (bar)")

    plt.subplot(122)
    plt.title("NH3 Abundance Models")

    #  NOMINAL ABUNDANCE
    # plt.plot(h2, p)
    # plt.plot(he, p)
    # plt.plot(ch4, p)
    plt.plot(nh3, p, label='NH3 $-$ Nominal')
    print nh3[p - 1 == np.min(np.abs(p - 1))]
    # plt.plot(h2o, p)
    # plt.legend(['h2', 'he', 'ch4', 'nh3', 'h2o'])
    tru_nh3 = np.array(nh3)



    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.e-8, 1.e-2])
    plt.ylim([1.e2, 1.e-2])
    plt.xlabel("Fractional Abundance")

    print '=='

    #  50 PCT
    # nh3 = np.copy(tru_nh3)
    # plim = np.where(p >= 0.5)
    # nh3[plim] *= 0.5
    # plt.plot(nh3, p, '--', label='_nolegend_')
    # plt.text(nh3[p == np.max(p)], 10*0.5, "50%")


    plt.show()


def run_sensitivity():
    f, tl, al, cl = mdl.examine_models(["0d5_60pct", "0d5_50pct", "0d5_45pct",
                                        "0d5_40pct", "0d5_30pct"])
    plt.figure(1)
    for i, a in enumerate(al):
        print "Abundance profile: {0}x at P>{1} bar".format(a, cl[i])
        plt.plot(f, tl[i], '.')
    plt.figure(2)
    tl = np.array(tl)
    # for i in range(tl.shape[1]):
    #     plt.plot()
    plt.show()


def run_contributions():
    color_wheel = paint_with_colors_of_wind()
    path_stub = "/home/ramsey/Documents/Research/pyplanet/Output/WGT/wgt_"
    txt = ".txt"
    alphas = [.3, .4, .45, .5, .6]

    plt.figure()
    # plt.subplot(121)

    for i in range(5):
        color = color_wheel.next()
        data = np.genfromtxt(path_stub + str(i) + txt, skip_header=1, usecols=[0, 2, 3])
        plt.plot(data[:, 1], data[:, 0], color=color, linestyle='--', label='_nolegend_')  # 27.5
        plt.plot(data[:, 2], data[:, 0], color=color, linestyle='-.', label=str(int(alphas[i]*100))+"%")  # 34.5
    plt.yscale('log')
    plt.ylim([2e0, 4e-1])
    plt.legend(loc='lower right')
    plt.title("Normalized Contribution Functions: 27.5 and 34.5 GHz")

    # plt.xscale('log')
    # plt.xlim([1e-10, 1e0])
    plt.plot(plt.xlim(), [1.8, 1.8], '--', color='black')
    plt.text(.03, 1.78, "P=1.8 bar")
    plt.plot(plt.xlim(), [0.46, 0.46], '--', color='black')
    plt.text(.03, 0.45, "P=0.46 bar")
    # plt.plot(plt.xlim(), [0.5, 0.5], '--', color='black')
    # plt.text(1e-6, 0.49, "P=0.5 bar")

    """
    CLOUDS vvv
    """
    """
    plt.subplot(122)

    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulclvla72x',
                         skip_header=1)

    layer = data[:, 0]
    t = data[:, 1]
    p = data[:, 2]
    soln = data[:, 3]
    h2o = data[:, 4]
    nh4sh = data[:, 5]
    nh3 = data[:, 6]

    def get_x_pos(x, y):
        # return x[x == np.max(x)]
        return gmean([x[x == np.max(x)], np.array([plt.xlim()[0]])])

    def get_y_pos(x, y):
        return np.sum(y[x != 0] * x[x != 0])/np.sum(x[x != 0])

    plt.plot(nh3, p, label='_nolegend_')
    # plt.text(get_x_pos(nh3, p), get_y_pos(nh3, p), "NH3")
    plt.plot(soln, p, label='_nolegend_')
    # plt.text(get_x_pos(soln, p), get_y_pos(soln, p), "SOLN")
    plt.plot(h2o, p, label='_nolegend_')
    # plt.text(get_x_pos(h2o, p), get_y_pos(h2o, p), "H2O")
    plt.plot(nh4sh, p, label='_nolegend_')
    # plt.text(get_x_pos(nh4sh, p), get_y_pos(nh4sh, p), "NH4SH")
    plt.xlabel("Cloud (?)")
    plt.ylabel("Pressure (bar)")
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([2e0, 4e-1])

    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlim([1.e-8, 1.e-2])
    # plt.ylim([1.e2, 1.e-2])
    # plt.xlabel("Fractional Abundance")
    plt.legend()

    """
    """
    CLOUDS ^^^
    """

    plt.tight_layout()
    plt.show()



run_atmosphere_graphs_2()
