"""
Quick run script for plotting T_B data
"""
import unpack as up
import plotting as pltt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.ticker import NullFormatter
import numpy as np
import datacontainer as dc
import models as mdl
import jackknife as jk
from itertools import cycle
from scipy.stats.mstats import gmean


def paint_with_colors_of_wind():
    return cycle(['cyan', 'blue', 'green', 'red',
                  'orange', 'maroon', 'black', 'purple'
                  ])


def run_normal():
    u = up.Unpack().prepare().adjust()
    wl = True
    p = pltt.Plotting(wl=wl)
    color_wheel = paint_with_colors_of_wind()
    model_group = "/home/ramsey/Documents/Research/pyplanet/Output/mgrid_physical_visualize/"
    for i, file_name in enumerate(up.generate_names(model_group)):
        t, r, d = mdl.examine_model(file_name)
        all_p = [t, r, d]
        count100 = 0
        weird_sum = t + r + d
        if (weird_sum != (79 + 108 + 42)) and (weird_sum != 300):
            continue
        label = "Top%d RH%d Deep%d" % (t, r, d) if weird_sum != 300 else "Nominal"
        p.add_model_plot(model_group + file_name, color=color_wheel.next(), line_sty='-',
                         label=label)
    # p.add_model_plot(mdl.PATH + stub1 + "/" + stub1 + ".dat", color=color_wheel.next(), line_sty='-',
    #                  label="1pct H, 35")
    # p.add_model_plot(mdl.PATH + stub2 + "/" + stub2 + ".dat", color=color_wheel.next(), line_sty='-',
    #                  label="10pct H, 35")
    # p.add_model_plot(mdl.PATH + stub3 + "/" + stub3 + ".dat", color=color_wheel.next(), line_sty='-',
    #                  label="10pct H, 45")
    # p.add_model_plot(mdl.PATH + stub4 + "/" + stub4 + ".dat", color=color_wheel.next(), line_sty='-',
    #                  label="1pct H, 45")
    p.add_my_data(u)
    p.add_all_other_data()
    u.write_points()
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


def run_large_emission_plot_wl():
    u = up.Unpack().prepare().adjust()
    wl = True
    p = pltt.Plotting(wl=wl)
    color_wheel = paint_with_colors_of_wind()

    def produce_name(top, rh, deep):
        model_name_1 = "/home/ramsey/Documents/Research/pyplanet/Output/final_models/t%drh" % top
        model_name_2 = "%dd%d.dat" % (rh, deep)
        return model_name_1 + model_name_2

    def produce_entire_model(top, rh, deep, color='k', linestyle='-', label=None, linewidth=1.):
        new_model = produce_name(top, rh, deep)
        p.add_model_plot(new_model, color=color, line_sty=linestyle, label=label, linewidth=linewidth)

    produce_entire_model(100, 100, 100, color='blue', linestyle='-', label="Nominal", linewidth=1.)
    produce_entire_model(22, 56, 42, color='red', linestyle='--', label="This work", linewidth=2.)
    produce_entire_model(15, 50, 39, color='black', linestyle='--', label="Low limits", linewidth=0.5)
    produce_entire_model(200, 64, 45, color='black', linestyle=':', label="High limits", linewidth=0.5)
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
    p.show_prep(location='lower left')
    p.show()


def run_large_emission_plot_f():
    u = up.Unpack().prepare().adjust()
    wl = False
    p = pltt.Plotting(wl=wl)
    fig = plt.figure(1)
    ax1 = plt.subplot(111)

    def produce_name(top, rh, deep):
        model_name_1 = "/home/ramsey/Documents/Research/pyplanet/Output/final_models/t%drh" % top
        model_name_2 = "%dd%d.dat" % (rh, deep)
        return model_name_1 + model_name_2

    def produce_entire_model(top, rh, deep, color='k', linestyle='-', label=None, linewidth=1., ax=None):
        new_model = produce_name(top, rh, deep)
        p.add_model_plot(new_model, color=color, line_sty=linestyle, label=label, linewidth=linewidth, ax=ax)

    produce_entire_model(100, 100, 100, color='blue', linestyle='-', label="Nominal", linewidth=1.)
    produce_entire_model(22, 56, 42, color='red', linestyle='--', label="This work", linewidth=2.)
    produce_entire_model(15, 50, 39, color='black', linestyle='--', label="Low limits", linewidth=0.5)
    produce_entire_model(200, 64, 45, color='black', linestyle=':', label="High limits", linewidth=0.5)
    produce_entire_model(100, 100, 42, color='cyan', linestyle='--', label="Adjusted Deep only", linewidth=1.)
    produce_entire_model(100, 56, 100, color='cyan', linestyle=':', label="Adjusted RH only", linewidth=1.)
    p.add_my_data(u)
    p.add_all_other_data()
    ax1.set_ylim([130, 180])
    ax1.set_ylabel("$T_{b}$ (K)")
    if wl:
        plt.title("Brightness Temperature (disk averaged) by Wavelength")
        plt.xlabel("$\lambda$ (cm)")
    else:
        plt.title("Brightness Temperature (disk averaged) by Frequency")
        plt.xlabel("$\\nu$ (GHz)")
        plt.xlim([100, 20])
    p.show_prep(location='upper right', ncol=3)
    inset_ax = inset_axes(ax1, width="40%", height="45%",
                          loc=3, borderpad=2)
    produce_entire_model(100, 100, 100, color='blue',
                         linestyle='-', label="Nominal", linewidth=1.,
                         ax=inset_ax)
    produce_entire_model(22, 56, 42, color='red',
                         linestyle='--', label="This work", linewidth=2.,
                         ax=inset_ax)
    produce_entire_model(15, 50, 39, color='black',
                         linestyle='--', label="Low limits", linewidth=0.5,
                         ax=inset_ax)
    produce_entire_model(200, 64, 45, color='black',
                         linestyle=':', label="High limits", linewidth=0.5,
                         ax=inset_ax)
    produce_entire_model(100, 100, 42, color='cyan', linestyle='--', label="Adjusted Deep only", linewidth=1.)
    produce_entire_model(100, 56, 100, color='cyan', linestyle=':', label="Adjusted RH only", linewidth=1.)
    p.add_my_data(u, ax=inset_ax)
    p.add_all_other_data(ax=inset_ax)
    inset_ax.set_xlim([36, 26])
    inset_ax.set_ylim([138, 155])
    inset_ax.yaxis.tick_right()
    plt.show()


def run_solar_vla():
    p = pltt.Plotting(wl=True)
    p.add_model_plot("paulsolar", color='orange', line_sty='-', linewidth=2.0)
    p.add_model_plot("paulvla72x", color='blue', line_sty='-', linewidth=2.0)
    p.show_prep()
    p.show()


def run_jackknife():
    default_u = up.Unpack().prepare().adjust()
    default_tb, default_tb_s, default_tb_o = default_u.get_temperatures()
    u = jk.JackknifeUnpack()
    plt.figure(10)
    counter = 1
    for t in jk.JACKKNIFE_TESTS:
        u.accept_test(jk.JACKKNIFE_TESTS[t][0]).adjust()
        tb, tbe_slope, tbe_offset = u.get_temperatures()
        tb -= default_tb
        mean_deviation = np.median(tb)
        spread_deviation = np.std(tb)
        plt.errorbar([counter], [mean_deviation], yerr=[spread_deviation],
                     fmt=jk.JACKKNIFE_TESTS[t][1],
                     color=jk.JACKKNIFE_TESTS[t][2],
                     label=t)
        counter += 1
    plt.title("Jackknife Tests")
    plt.ylabel("Variation from $T_{b}$ (K)")
    plt.xlabel("Test #")
    plt.legend()
    plt.show()
    # b = jk.BootstrapAverageContainer(f[:f.size / 2], u)
    # b.cycle(data.tb)
    # plt.plot(b.f, b.t, ',', color='k', label="Bootstrap")
    # plt.ylim([0 - 2 * max_error, 0 + 2 * max_error])


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
        plt.errorbar(u.get_dates(), ch.flux, yerr=ch.error, fmt='.')
        print "F %.3f: %d" % (ch.frequency_ghz, len(np.where(np.isfinite(ch.flux))[0]))
    plt.xlim(date_min - l_padding * date_range, np.max(u.dates_copss_array) + l_padding * date_range)
    plt.subplot(122)
    plt.title("Raw Data - After Normalization")
    plt.xlabel("Date (MJD)")
    # plt.ylabel("Flux density (jy)")
    u.distance_adj()
    for ch in u.channel_obj_list:
        plt.errorbar(u.get_dates(), ch.flux, yerr=ch.error, fmt='.')
        plt.text(date_max + 3,
                 np.nanmean(ch.flux),
                 "%.3f GHz" % ch.frequency_ghz, fontsize=8)
    plt.xlim(date_min - l_padding * date_range, np.max(u.dates_copss_array) + r_padding * date_range)
    plt.show()


def run_raw_airmass():
    """
    Figure-ready, raw data before and after airmass correction
    """
    u = up.Unpack().prepare()
    plt.figure()
    # plt.subplot(121)
    plt.title("Raw Data - Before Normalization")
    plt.xlabel("Date (MJD)")
    plt.ylabel("Flux density (jy)")
    date_max = np.max(u.dates_copss_array)
    date_min = np.min(u.dates_copss_array)
    date_range = date_max - date_min
    l_padding = 0.05
    r_padding = 0.2
    u.distance_adj()
    stds = {}
    for ch in u.channel_obj_list:
        plt.errorbar(u.get_j_el(), ch.flux, yerr=ch.error, fmt='x')
        # plt.errorbar(u.dates_copss_array, ch.flux, yerr=ch.error, fmt='x')
        stds[str(ch.frequency_ghz)] = np.nanstd(ch.flux)
        flux = ch.flux[np.where(ch.flux < np.inf)]
        j_el = u.get_j_el()[np.where(ch.flux < np.inf)]
        e_fit = np.polyfit(j_el, flux, deg=1)
        e_fit_line = mdl.polynomial(e_fit, j_el)
        plt.plot(j_el, e_fit_line, '-', )
    # plt.xlim(date_min - l_padding * date_range, np.max(u.dates_copss_array) + l_padding * date_range)
    # plt.subplot(122)
    plt.title("Raw Data - After Normalization")
    plt.xlabel("Date (MJD)")
    plt.ylabel("Flux density (jy)")
    u.alt_adj()
    for ch in u.channel_obj_list:
        plt.errorbar(u.get_j_el(), ch.flux, yerr=ch.error, fmt='.')
        # plt.errorbar(u.dates_copss_array, ch.flux, yerr=ch.error, fmt='.')
        stds[str(ch.frequency_ghz)] = 100 * (stds[str(ch.frequency_ghz)] - np.nanstd(ch.flux)) / stds[
            str(ch.frequency_ghz)]
        flux = ch.flux[np.where(ch.flux < np.inf)]
        j_el = u.get_j_el()[np.where(ch.flux < np.inf)]
        e_fit = np.polyfit(j_el, flux, deg=1)
        e_fit_line = mdl.polynomial(e_fit, j_el)
        plt.plot(j_el, e_fit_line, '--', )
        # plt.text(date_max + 3,
        #          np.nanmean(ch.flux),
        #          "%.3f GHz" % ch.frequency_ghz, fontsize=8)
    # plt.xlim(date_min - l_padding * date_range, np.max(u.dates_copss_array) + r_padding * date_range)
    pairs = [(k, stds[k]) for k in stds]
    f_list, std_list = zip(*pairs)
    f_list = np.array(f_list)
    std_list = np.array(std_list)
    ordered = np.argsort(f_list)
    f_list = f_list[ordered]
    std_list = std_list[ordered]
    for f, s in zip(f_list, std_list):
        print "F", f, "delta std", s
    plt.show()


def run_get_relative_uncertainty():
    """
    Trying to get relative uncertainty
    """
    u = up.Unpack().prepare()
    u.distance_adj().alt_adj().prelim_adjust().average_adj()
    u.error_adj()
    plt.figure()
    plt.title("Relative Uncertainty")
    flux_grid = np.array([]).reshape((0, 37))
    flux_error_grid = np.array([]).reshape((0, 37))
    flux_avg_vector = np.array([])
    frequency_vector = np.array([])
    frq_axis = 1
    day_axis = 0
    for ch in u.channel_obj_list:
        flux_grid = np.concatenate([flux_grid, np.array([ch.flux])])
        flux_error_grid = np.concatenate([flux_error_grid, np.array([ch.error])])
        flux_avg_vector = np.append(flux_avg_vector, ch.flux_avg)
        frequency_vector = np.append(frequency_vector, ch.frequency_ghz)
    flux_avg_vector = np.transpose(np.array([flux_avg_vector]))
    frequency_vector = np.transpose(np.array([frequency_vector]))
    fractional_error_by_channel = np.nanmean(flux_error_grid / flux_grid, axis=frq_axis) * flux_avg_vector.transpose()
    fractional_error_by_channel = np.transpose(fractional_error_by_channel)
    flux_grid /= flux_avg_vector
    # weighting was a bad idea
    # flux_grid /= flux_error_grid**2.
    # weighting was a bad idea
    # flux_grid /= np.reciprocal(np.nanmin(flux_error_grid, axis=day_axis)).reshape([1, 37])**2.
    average_by_day = np.nanmean(flux_grid, axis=day_axis)
    spread_by_day = np.nanstd(flux_grid, axis=day_axis)

    # daily spread!
    # plt.plot(u.get_dates(), flux_grid.transpose(), '.')
    # plt.errorbar(u.get_dates(), average_by_day, yerr=spread_by_day, fmt='.', color='black', zorder=1000)
    # plt.xlabel("Obs Date (MJD)")
    # plt.ylabel("Normalized Measured Flux")

    flux_grid -= average_by_day
    flux_grid *= flux_avg_vector
    # flux_grid /= average_by_day
    # flux_grid -= 1

    # daily spread centered at 0
    # plt.plot(u.get_dates(), flux_grid.transpose(), '.')
    # plt.errorbar(u.get_dates(), np.zeros(u.get_dates().size), yerr=spread_by_day, fmt='.', color='black', zorder=1000)
    # plt.xlabel("Obs Date (MJD)")
    # plt.ylabel("Normalized Measured Flux")

    spread_daily_deviation_by_channel = np.nanstd(flux_grid, axis=frq_axis)
    # flux_grid /= spread_by_day  # put it in terms of stds
    average_daily_deviation_by_channel = np.nanmean(flux_grid, axis=frq_axis)

    # Deviation & spread!
    # plt.plot(frequency_vector, flux_grid)  # THIS WARRANTS A QUESTION
    # plt.errorbar(frequency_vector, average_daily_deviation_by_channel,
    #              yerr=spread_daily_deviation_by_channel, fmt='o')
    # plt.xlabel("Frequency (GHz)")
    # plt.ylabel("Mean Deviation from Normalized Daily Average (jy)")

    spread_error_by_channel = spread_daily_deviation_by_channel.reshape((15, 1))
    deviation_error_by_channel = average_daily_deviation_by_channel.reshape((15, 1))

    # check order of relative error contributions
    plt.plot(frequency_vector, np.sqrt(fractional_error_by_channel ** 2.), 'x', label="Fractional")
    plt.plot(frequency_vector, np.sqrt(deviation_error_by_channel ** 2.), '^', label="Mean Deviation")
    plt.plot(frequency_vector, np.sqrt(spread_error_by_channel ** 2.), 'D', label="Spread Deviation")
    plt.yscale('log')

    # FINAL VALUE !!!!!!!!!!!!!!!!!!!!!!!!!!
    relative_uncertainty = fractional_error_by_channel ** 2 + deviation_error_by_channel ** 2 + spread_error_by_channel ** 2
    plt.plot(frequency_vector, np.sqrt(relative_uncertainty), '*', label="Relative")
    print (np.sqrt(relative_uncertainty) * 100 / frequency_vector).reshape(15)
    print np.sqrt(relative_uncertainty).reshape(15)

    # something else interesting kinda sorta
    # flux_grid_errors = np.nanstd(flux_grid, axis=a) / np.nansum(np.reciprocal(flux_error_grid**2.), axis=a)
    # flux_grid_averages = np.nansum(flux_grid, axis=a) / np.nansum(np.reciprocal(flux_error_grid**2.), axis=a)
    # plt.plot(u.get_frequencies(), flux_grid_errors, '.')
    # plt.errorbar(u.get_frequencies(), flux_grid_averages, yerr=flux_grid_errors, fmt='.')
    # plt.plot(frequency_vector, flux_avg_vector/flux_avg_vector, 'o')

    # trying to look at channel spread
    # plt.plot(frequency_vector, flux_grid, 'x')
    # plt.plot(frequency_vector, flux_error_grid, 'x')

    # WEIRD PARABOLA was from error weighting......

    # plt.xscale('log')
    plt.legend()
    plt.show()


def run_model_grid_2d_first_iteration():
    unpacker = up.Unpack().prepare().adjust()
    tb, tbe_slope, tbe_offset = unpacker.get_temperatures()
    frequencies = unpacker.get_frequencies()
    data_model_info = (frequencies, tb, tbe_slope, tbe_offset)
    all_others = True
    #  param 1 is relative humidity
    #  param 2 is deep abundance
    (s1, e1, step1) = (.9, 1.12, .02)
    (s2, e2, step2) = (.24, .45, .2)
    params1, chi_sq_values1 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_vary_rh_deep_1/",
                                                         all_others=all_others)
    (s1, e1, step1) = (.92, 1.12, .04)
    (s2, e2, step2) = (.3, .45, .02)
    params2, chi_sq_values2 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_vary_rh_deep_2/",
                                                         all_others=all_others)
    (s1, e1, step1) = (1.05, 1.11, .01)
    (s2, e2, step2) = (.35, .48, .01)
    params3, chi_sq_values3 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_vary_rh_deep_3/",
                                                         all_others=all_others)
    params = np.concatenate([params1, params2, params3])
    chi_sq_values = np.concatenate([chi_sq_values1, chi_sq_values2, chi_sq_values3])
    humidity = np.array([x[0] for x in params]) * 100
    deep_abundance = np.array([x[1] for x in params]) * 100
    lowest_chsq_array = np.array(chi_sq_values)
    lowest_val = np.where(lowest_chsq_array == np.min(lowest_chsq_array))
    lowest_param_array = [params[i] for i in lowest_val[0]]
    print "Chi Sq Min = %.3f" % np.min(lowest_chsq_array)
    print "Lowest Chi Squared Values:"
    for (p1, p2) in lowest_param_array:
        print "RH: %.2f, Deep:%.2f" % (round(p1, 3) * 100, round(p2, 3) * 100)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(
        deep_abundance, humidity, zs=chi_sq_values,
        color='black', marker='o', linestyle='None'
    )
    # ax.set_xlim([35, 40])
    # ax.set_ylim([108, 112])
    # ax.set_zlim([0, 2])  # all_others = True
    plt.xlabel("Deep Abundance %")
    plt.ylabel("Rel Humidity %")
    ax.set_zlabel("Reduced $\chi^{2}$")
    plt.show()


def run_model_grid_2d_second_iteration():
    unpacker = up.Unpack().prepare().adjust()
    tb, tbe_slope, tbe_offset = unpacker.get_temperatures()
    frequencies = unpacker.get_frequencies()
    data_model_info = (frequencies, tb, tbe_slope, tbe_offset)
    all_others = True
    #  param 1 is relative humidity
    #  param 2 is deep abundance

    # ALL
    (s1, e1, step1) = (1.0, 1.15, .02)
    (s2, e2, step2) = (.35, .51, .02)
    params1, chi_sq_values1 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_vary_rh_deep_top76_1/",
                                                         all_others=all_others)
    (s1, e1, step1) = (1.07, 1.14, .02)
    (s2, e2, step2) = (.36, .44, .02)
    params2, chi_sq_values2 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_vary_rh_deep_top76_2/",
                                                         all_others=all_others)
    (s1, e1, step1) = (1.06, 1.12, .02)
    (s2, e2, step2) = (.36, .44, .02)
    params3, chi_sq_values3 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_vary_rh_deep_top76_3/",
                                                         all_others=all_others)

    params = np.concatenate([params1, params2, params3])
    chi_sq_values = np.concatenate([chi_sq_values1, chi_sq_values2, chi_sq_values3])
    humidity = np.array([x[0] for x in params]) * 100
    deep_abundance = np.array([x[1] for x in params]) * 100
    lowest_chsq_array = np.array(chi_sq_values)
    lowest_val = np.where(lowest_chsq_array == np.min(lowest_chsq_array))
    lowest_param_array = [params[i] for i in lowest_val[0]]
    print "Chi Sq Min = %.3f" % np.min(lowest_chsq_array)
    print "Lowest Chi Squared Values:"
    for (p1, p2) in lowest_param_array:
        print "RH: %.2f, Deep:%.2f" % (round(p1, 3) * 100, round(p2, 3) * 100)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(
        deep_abundance, humidity, zs=chi_sq_values,
        color='black', marker='o', linestyle='None'
    )
    # ax.set_xlim([41, 43])
    # ax.set_ylim([107, 109])
    # ax.set_zlim([0, 1])  # all_others = True
    plt.xlabel("Deep Abundance %")
    plt.ylabel("Rel Humidity %")
    ax.set_zlabel("Reduced $\chi^{2}$")
    plt.show()


def run_model_grid_2d_third_iteration():  # This goes back to multiplicative RH s.t. constant thru clouds
    """
    Fixed for RH multiplicative 7/9/17
    :return:
    """
    unpacker = up.Unpack().prepare().adjust()
    tb, tbe_slope, tbe_offset = unpacker.get_temperatures()
    frequencies = unpacker.get_frequencies()
    data_model_info = (frequencies, tb, tbe_slope, tbe_offset)
    all_others = True
    #  param 1 is relative humidity
    #  param 2 is deep abundance

    # ALL
    (s1, e1, step1) = (.45, .65, .02)
    (s2, e2, step2) = (.38, .48, .02)
    params1, chi_sq_values1 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_multrh_vary_rh_deep_top76_1/",
                                                         all_others=all_others)
    (s1, e1, step1) = (.54, .561, 0.005)
    (s2, e2, step2) = (.41, .431, 0.005)
    params2, chi_sq_values2 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
                                                         data_model_info,
                                                         "mgrid_multrh_vary_rh_deep_top76_2/",
                                                         all_others=all_others)
    # (s1, e1, step1) = (1.06, 1.12, .02)
    # (s2, e2, step2) = (.36, .44, .02)
    # params3, chi_sq_values3 = mdl.create_parameter_space((s1, e1, step1), (s2, e2, step2),
    #                                                      data_model_info,
    #                                                      "mgrid_vary_rh_deep_top76_3/",
    #                                                      all_others=all_others)
    #
    params = np.concatenate([params1, params2])
    chi_sq_values = np.concatenate([chi_sq_values1, chi_sq_values2])
    humidity = np.array([x[0] for x in params]) * 100
    deep_abundance = np.array([x[1] for x in params]) * 100
    lowest_chsq_array = np.array(chi_sq_values)
    lowest_val = np.where(lowest_chsq_array == np.min(lowest_chsq_array))
    lowest_param_array = [params[i] for i in lowest_val[0]]
    print "Chi Sq Min = %.3f" % np.min(lowest_chsq_array)
    print "Lowest Chi Squared Values:"
    for (p1, p2) in lowest_param_array:
        print "RH: %.2f, Deep:%.2f" % (round(p1, 3) * 100, round(p2, 3) * 100)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(
        deep_abundance, humidity, zs=chi_sq_values,
        color='black', marker='o', linestyle='None'
    )
    # ax.set_xlim([41, 43])
    # ax.set_ylim([107, 109])
    # ax.set_zlim([0, 1])  # all_others = True
    plt.xlabel("Deep Abundance %")
    plt.ylabel("Rel Humidity %")
    ax.set_zlabel("Reduced $\chi^{2}$")
    plt.show()


def run_model_grid_1d_first_iteration():
    """
    Runs ONE parameter.
    Remodeled for top abundance 6/22/17
    :return:
    """
    deg = 8  # To get a good interpolation of model fit
    unpacker = up.Unpack().prepare().adjust()
    tb, tbe_slope, tbe_offset = unpacker.get_temperatures()
    frequencies = unpacker.get_frequencies()
    data_model_info = (frequencies, tb, tbe_slope, tbe_offset)
    color_wheel = paint_with_colors_of_wind()
    markers = ['o', 's', 'v', '^', '*']
    labels = ['relative only', 'absolute + WMAP + Gibson']
    x_lim = [-.1, 2.5]
    all_others = True
    # parameter is top abundance
    # ALL
    params1, chi_squared_values1 = mdl.create_parameter_space_1d((0.1, 1.5, 0.1),
                                                                 data_model_info, "mgrid_vary_top_rh108deep42_1/",
                                                                 all_others=all_others)  # RH108 Deep42, per ALL
    params2, chi_squared_values2 = mdl.create_parameter_space_1d(np.concatenate([np.array([.001, .01, .05]),
                                                                                 np.arange(1.5, 2.5, .3)]),
                                                                 data_model_info, "mgrid_vary_top_rh108deep42_2/",
                                                                 all_others=all_others)  # RH108 Deep42, per ALL
    params3, chi_squared_values3 = mdl.create_parameter_space_1d(10. ** np.arange(-6, -3, 1),
                                                                 data_model_info, "mgrid_vary_top_rh108deep42_3/",
                                                                 all_others=all_others)  # RH108 Deep42, per ALL
    params4, chi_squared_values4 = mdl.create_parameter_space_1d(np.arange(2.7, 3.7, .3),
                                                                 data_model_info, "mgrid_vary_top_rh108deep42_4/",
                                                                 all_others=all_others)  # RH108 Deep42, per ALL
    params = np.concatenate([params1, params2, params3, params4])
    chi_squared_values = np.concatenate([chi_squared_values1, chi_squared_values2,
                                         chi_squared_values3, chi_squared_values4])

    # SLOPE
    # params1, chi_squared_values1 = mdl.create_parameter_space_1d(np.concatenate([np.array([.001, .01, .05]),
    #                                                                              np.arange(0.1, 1.5, 0.1),
    #                                                                              np.arange(1.5, 2.5, .3)]),
    #                                                              data_model_info, "mgrid_vary_top_rh110deep36_1/",
    #                                                              all_others=all_others)  # RH110 Deep36, per SLOPE
    # params2, chi_squared_values2 = mdl.create_parameter_space_1d((.1, .3, .01),
    #                                                              data_model_info, "mgrid_vary_top_rh110deep36_2/",
    #                                                              all_others=all_others)  # RH110 Deep36, per SLOPE
    # params = np.concatenate([params1, params2])
    # chi_squared_values = np.concatenate([chi_squared_values1, chi_squared_values2])

    clr = color_wheel.next()
    top_abundance = np.array(params)
    plt.plot(params, chi_squared_values,
             markers[all_others], color=clr,
             label=labels[all_others])

    fit = np.polyfit(top_abundance, chi_squared_values, deg=deg)
    x_range = np.arange(x_lim[0], x_lim[1], 0.01)
    output = mdl.polynomial(fit, x_range)
    plt.plot(x_range, output, '--', color=clr)
    where_min = np.argmin(output)
    top_ab_minimum = x_range[where_min]
    fit_errors = [0, 0]
    # noinspection PyTypeChecker
    where_lo_error = np.where(
        np.min((output[where_min] * 2. - output[:where_min]) ** 2) == (
            output[where_min] * 2. - output[:where_min]) ** 2)
    fit_errors[0] = top_ab_minimum - x_range[:where_min][where_lo_error]
    # noinspection PyTypeChecker
    where_hi_error = np.where(
        np.min((output[where_min:] - output[where_min] * 2.) ** 2) == (
            output[where_min:] - output[where_min] * 2.) ** 2)
    fit_errors[1] = x_range[where_min:][where_hi_error] - top_ab_minimum
    # plt.errorbar([top_ab_minimum], [output[where_min] * 1.01],
    #              xerr=[fit_errors[0], fit_errors[1]],
    #              fmt=markers[all_others],
    #              color=clr, label=labels[all_others])
    print "MIN", top_ab_minimum
    # plt.xscale('log')
    plt.xlabel("Top Abundance Multiplier")
    plt.ylabel("Reduced $\chi^{2}$")
    plt.title("Model comparison to CARMA and WMAP data")
    plt.legend(loc='upper right')
    plt.show()


def run_model_grid_1d_second_iteration():  # This is after RH changed back to multiplicative
    """
    Runs ONE parameter.
    Remodeled for top abundance 6/22/17
    Fixed for RH multiplicative 7/9/17
    :return:
    """
    deg = 5  # To get a good interpolation of model fit
    unpacker = up.Unpack().prepare().adjust()
    tb, tbe_slope, tbe_offset = unpacker.get_temperatures()
    frequencies = unpacker.get_frequencies()
    data_model_info = (frequencies, tb, tbe_slope, tbe_offset)
    color_wheel = paint_with_colors_of_wind()
    markers = ['o', 's', 'v', '^', '*']
    labels = ['relative only', 'absolute + WMAP + Gibson']
    x_lim = [.01, 2.4]
    all_others = True
    # parameter is top abundance
    # ALL
    params1, chi_squared_values1 = mdl.create_parameter_space_1d((.7, 1.11, .1),
                                                                 data_model_info, "mgrid_multrh_vary_top_rh55deep42_1/",
                                                                 all_others=all_others)  # RH55 Deep42, per ALL
    params2, chi_squared_values2 = mdl.create_parameter_space_1d(np.arange(0.1, 0.61, .2),
                                                                 data_model_info, "mgrid_multrh_vary_top_rh55deep42_2/",
                                                                 all_others=all_others)  # RH55 Deep42, per ALL
    params3, chi_squared_values3 = mdl.create_parameter_space_1d(10. ** np.arange(-6, -2, 1),
                                                                 data_model_info, "mgrid_multrh_vary_top_rh55deep42_3/",
                                                                 all_others=all_others)  # RH108 Deep42, per ALL
    params4, chi_squared_values4 = mdl.create_parameter_space_1d(np.arange(1.2, 2.5, .4),
                                                                 data_model_info, "mgrid_multrh_vary_top_rh55deep42_4/",
                                                                 all_others=all_others)  # RH108 Deep42, per ALL
    params = np.concatenate([params1, params2, params3, params4])
    chi_squared_values = np.concatenate([chi_squared_values1, chi_squared_values2,
                                         chi_squared_values3, chi_squared_values4])

    clr = color_wheel.next()
    top_abundance = np.array(params)
    plt.plot(params, chi_squared_values,
             markers[all_others], color=clr,
             label=labels[all_others])

    fit = np.polyfit(top_abundance, chi_squared_values, deg=deg)
    x_range = np.arange(x_lim[0], x_lim[1], 0.01)
    output = mdl.polynomial(fit, x_range)
    plt.plot(x_range, output, '--', color=clr)
    where_min = np.argmin(output)
    top_ab_minimum = x_range[where_min]
    fit_errors = [0, 0]
    # noinspection PyTypeChecker
    where_lo_error = np.where(
        np.min((output[where_min] * 2. - output[:where_min]) ** 2) == (
            output[where_min] * 2. - output[:where_min]) ** 2)
    fit_errors[0] = top_ab_minimum - x_range[:where_min][where_lo_error]
    # noinspection PyTypeChecker
    where_hi_error = np.where(
        np.min((output[where_min:] - output[where_min] * 2.) ** 2) == (
            output[where_min:] - output[where_min] * 2.) ** 2)
    fit_errors[1] = x_range[where_min:][where_hi_error] - top_ab_minimum
    plt.errorbar([top_ab_minimum], [output[where_min] * 2.],
                 xerr=[fit_errors[0], fit_errors[1]],
                 fmt=markers[all_others],
                 color=clr, label=labels[all_others])
    print "MIN %.1f : [%.1f, %.1f]" % (top_ab_minimum * 100, 100 * (top_ab_minimum - fit_errors[0]),
                                       100 * (top_ab_minimum + fit_errors[1]))
    # plt.xscale('log')
    plt.xlabel("Top Abundance Multiplier")
    plt.ylabel("Reduced $\chi^{2}$")
    plt.title("Model comparison to CARMA and WMAP data")
    plt.legend(loc='upper right')
    plt.show()


def run_model_grid_1d_errorbars_rh():
    deg = 8  # To get a good interpolation of model fit
    unpacker = up.Unpack().prepare().adjust()
    tb, tbe_slope, tbe_offset = unpacker.get_temperatures()
    frequencies = unpacker.get_frequencies()
    data_model_info = (frequencies, tb, tbe_slope, tbe_offset)
    color_wheel = paint_with_colors_of_wind()
    markers = ['o', 's', 'v', '^', '*']
    labels = ['relative only', 'absolute + WMAP + Gibson']
    x_lim = [.4, .67]
    all_others = True
    params1, chi_squared_values1 = mdl.create_parameter_space_1d((.4, .67, 0.05),
                                                                 data_model_info, "mgrid_multrh_vary_rh_top22deep42_1/",
                                                                 all_others=all_others)
    params2, chi_squared_values2 = mdl.create_parameter_space_1d((.415, .67, 0.05),
                                                                 data_model_info, "mgrid_multrh_vary_rh_top22deep42_2/",
                                                                 all_others=all_others)
    params = np.concatenate([params1, params2])
    chi_squared_values = np.concatenate([chi_squared_values1, chi_squared_values2])
    clr = color_wheel.next()
    top_abundance = np.array(params)
    plt.plot(params, chi_squared_values,
             markers[all_others], color=clr,
             label=labels[all_others])

    fit = np.polyfit(top_abundance, chi_squared_values, deg=deg)
    x_range = np.arange(x_lim[0], x_lim[1], 0.005)
    output = mdl.polynomial(fit, x_range)
    plt.plot(x_range, output, '--', color=clr)
    where_min = np.argmin(output)
    top_ab_minimum = x_range[where_min]
    fit_errors = [0, 0]
    # noinspection PyTypeChecker
    where_lo_error = np.where(
        np.min((output[where_min] * 2. - output[:where_min]) ** 2) == (
            output[where_min] * 2. - output[:where_min]) ** 2)
    fit_errors[0] = top_ab_minimum - x_range[:where_min][where_lo_error]
    # noinspection PyTypeChecker
    where_hi_error = np.where(
        np.min((output[where_min:] - output[where_min] * 2.) ** 2) == (
            output[where_min:] - output[where_min] * 2.) ** 2)
    fit_errors[1] = x_range[where_min:][where_hi_error] - top_ab_minimum
    plt.errorbar([top_ab_minimum], [output[where_min] * 3],
                 xerr=[fit_errors[0], fit_errors[1]],
                 fmt=markers[all_others],
                 color=clr, label=labels[all_others])
    print "MIN %.1f : [%.1f, %.1f]" % (top_ab_minimum * 100, 100 * (top_ab_minimum - fit_errors[0]),
                                       100 * (top_ab_minimum + fit_errors[1]))
    plt.xlabel("Relative Humidity Multiplier")
    plt.ylabel("Reduced $\chi^{2}$")
    plt.title("Model comparison to CARMA and WMAP data")
    plt.legend(loc='upper right')
    plt.show()


def run_model_grid_1d_errorbars_deep():
    deg = 8  # To get a good interpolation of model fit
    unpacker = up.Unpack().prepare().adjust()
    tb, tbe_slope, tbe_offset = unpacker.get_temperatures()
    frequencies = unpacker.get_frequencies()
    data_model_info = (frequencies, tb, tbe_slope, tbe_offset)
    color_wheel = paint_with_colors_of_wind()
    markers = ['o', 's', 'v', '^', '*']
    labels = ['relative only', 'absolute + WMAP + Gibson']
    x_lim = [.25, .6]
    all_others = True
    params1, chi_squared_values1 = mdl.create_parameter_space_1d((.25, .6, 0.02),
                                                                 data_model_info, "mgrid_multrh_vary_deep_top22rh55_1/",
                                                                 all_others=all_others)
    params = params1
    chi_squared_values = chi_squared_values1
    clr = color_wheel.next()
    top_abundance = np.array(params)
    plt.plot(params, chi_squared_values,
             markers[all_others], color=clr,
             label=labels[all_others])

    fit = np.polyfit(top_abundance, chi_squared_values, deg=deg)
    x_range = np.arange(x_lim[0], x_lim[1], 0.005)
    output = mdl.polynomial(fit, x_range)
    plt.plot(x_range, output, '--', color=clr)
    where_min = np.argmin(output)
    top_ab_minimum = x_range[where_min]
    fit_errors = [0, 0]
    # noinspection PyTypeChecker
    where_lo_error = np.where(
        np.min((output[where_min] * 2. - output[:where_min]) ** 2) == (
            output[where_min] * 2. - output[:where_min]) ** 2)
    fit_errors[0] = top_ab_minimum - x_range[:where_min][where_lo_error]
    # noinspection PyTypeChecker
    where_hi_error = np.where(
        np.min((output[where_min:] - output[where_min] * 2.) ** 2) == (
            output[where_min:] - output[where_min] * 2.) ** 2)
    fit_errors[1] = x_range[where_min:][where_hi_error] - top_ab_minimum
    plt.errorbar([top_ab_minimum], [output[where_min] * 20],
                 xerr=[fit_errors[0], fit_errors[1]],
                 fmt=markers[all_others],
                 color=clr, label=labels[all_others])
    print "MIN %.1f : [%.1f, %.1f]" % (top_ab_minimum * 100, 100 * (top_ab_minimum - fit_errors[0]),
                                       100 * (top_ab_minimum + fit_errors[1]))
    plt.xlabel("Deep Abundance Multiplier")
    plt.ylabel("Reduced $\chi^{2}$")
    plt.title("Model comparison to CARMA and WMAP data")
    plt.legend(loc='upper right')
    plt.show()


def run_air_mass_investigation():
    u = up.Unpack().prepare()
    u.distance_adj()
    # u.alt_adj()
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
    plt.xlim([24, 79])
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
    s_fit_line = mdl.polynomial(s_fit, freqs)
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

    plt.show(11)


def run_atmosphere_graphs():
    file_stub = "/home/ramsey/Documents/Research/pyplanet/Jupiter/for_ramsey/"

    # TEMPERATURE-PRESSURE
    plt.figure(0)
    # TP Graph
    plt.subplot(121)
    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulvla72x',
                         skip_header=1)
    plt.plot(data[:, 1], data[:, 2], color='blue', linestyle='-', label="VLA 72x")
    # plt.subplot(122)
    # plt.plot(data[:, 2] * 1.342E7 * np.exp(-3753.6 / data[:, 1]),
    #          data[:, 2],  # Clausius Clapeyron
    #          color='black', linewidth=1
    #          )
    # plt.subplot(121)
    # data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulvla72',
    #                      skip_header=1)
    # plt.plot(data[:, 1], data[:, 2], color='green', linestyle='-.', label="VLA 72")

    # plt.subplot(122)
    # plt.plot(data[:, 2] * 1.342E7 * np.exp(-3753.6 / data[:, 1]),
    #          data[:, 2],  # Clausius Clapeyron
    #          color='black', linewidth=1
    #          )
    # plt.subplot(121)

    # plot logistics
    plt.title("Temperature-Pressure Profile")
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1.e1, 1.e-1])
    plt.xlim(([1.e2, 1.e3]))
    plt.xlabel("T (K)")
    plt.ylabel("P (bar)")
    # plt.legend()

    """
    Contributions
    *EDIT: we need to run new ones RESOLVED
    **RESOLVED: we ran new ones; they look GOOD
    """
    txt = ".txt"
    linestys = ['-.', '--', '-', ':']
    wgt_freqs = [22.5, 27.5, 34.5, 90.5]
    wgt_nom = np.genfromtxt(file_stub + "wgt_nominal" + txt, skip_header=1, usecols=[0, 2, 3, 4, 5])
    wgt_t22rd55d42 = np.genfromtxt(file_stub + "wgt_t22rh55d42" + txt, skip_header=1, usecols=[0, 2, 3, 4, 5])
    anchor_x = 3
    reach_x = .5
    for i in range(4):
        i_text = '%.1f GHz' % wgt_freqs[i]
        plt.plot(10 ** (anchor_x - reach_x * wgt_nom[:, i + 1]), wgt_nom[:, 0],
                 color='blue', linestyle=linestys[i], label='_nolegend_',
                 linewidth=0.5)
        plt.text(10 ** (anchor_x - 1.2 * reach_x * wgt_nom[:, i + 1][wgt_nom[:, i + 1] == np.max(wgt_nom[:, i + 1])]),
                 wgt_nom[:, 0][wgt_nom[:, i + 1] == np.max(wgt_nom[:, i + 1])],
                 i_text, fontsize=8
                 )
        plt.plot(10 ** (anchor_x - reach_x * wgt_t22rd55d42[:, i + 1]), wgt_t22rd55d42[:, 0],
                 color='red', linestyle=linestys[i], label='_nolegend_',
                 linewidth=0.5)

    # CONSTITUENT ABUNDANCES
    plt.subplot(122)
    plt.title("NH3 Abundance Models")

    def plot_profile(name, color='red', linestyle='-', label='_nolegend_', linewidth=1.):
        nh3_abundance = np.loadtxt(file_stub + name)
        plt.plot(nh3_abundance[1, :], nh3_abundance[0, :],
                 color=color, linestyle=linestyle, label=label, linewidth=linewidth)

    plot_profile("nh3_abundance_nominal.dat", color='blue', linestyle='-',
                 label="NH3$-\ P\sim 0.1: 1.2\\times 10^{-7} \ -\ 100\% RH \ -\ P\sim 8: 5.72\\times 10^{-4}$"
                       " (Nominal)",
                 linewidth=1)
    plot_profile("nh3_abundance_t22rh56d42.dat", color='red', linestyle='--',
                 label="NH3$-\ P\sim 0.1: 2.8\\times 10^{-8} \ -\ 55\% RH \ -\ P\sim 8: 2.40\\times 10^{-4}$"
                       " (This work)",
                 linewidth=2)
    plot_profile("nh3_abundance_t15rh50d39.dat", color='black', linestyle='--',
                 label="NH3$-\ P\sim 0.1: 1.9\\times 10^{-8} \ -\ 50\% RH \ -\ P\sim 8: 2.26\\times 10^{-4}$"
                       " (Low limits)",
                 linewidth=0.5)
    plot_profile("nh3_abundance_t200rh64d45.dat", color='black', linestyle=':',
                 label="NH3$-\ P\sim 0.1: 2.4\\times 10^{-7} \ -\ 64\% RH \ -\ P\sim 8: 2.57\\times 10^{-4}$"
                       " (High limits)",
                 linewidth=0.5)

    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulvla72',
                         skip_header=1)
    t = data[:, 1]
    p = data[:, 2]
    nh3 = data[:, 6]
    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulSolar_asplund_ws',
                         skip_header=1)
    ps_deep = data[:, 2]
    nh3s_deep = data[:, 6][ps_deep > 0.5]
    ps_deep = ps_deep[ps_deep > 0.5]
    ps = np.concatenate([p[p < 0.5], ps_deep])
    nh3s = np.concatenate([nh3[p < 0.5], nh3s_deep])
    plt.plot(nh3s, ps, color="orange", linestyle="-", label='Solar', linewidth=1)
    # plot logistics
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.e-8, 1.e-2])
    plt.ylim([1.e1, 1.e-1])
    plt.xlabel("Fractional Abundance")

    # CLOUDS
    data = np.genfromtxt('/home/ramsey/Documents/Research/pyplanet/Jupiter/jupiter.paulclvla72x',
                         skip_header=1)
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
        return np.sum(y[x != 0] * x[x != 0]) / np.sum(x[x != 0])

    plt.plot(nh3, p, label='_nolegend_', linewidth=0.5)
    plt.text(get_x_pos(nh3, p), get_y_pos(nh3, p), "NH3")
    plt.plot(soln, p, label='_nolegend_', linewidth=0.5)
    plt.text(get_x_pos(soln, p), get_y_pos(soln, p), "SOLN")
    plt.plot(h2o, p, label='_nolegend_', linewidth=0.5)
    plt.text(get_x_pos(h2o, p), get_y_pos(h2o, p), "H2O")
    plt.plot(nh4sh, p, label='_nolegend_', linewidth=0.5)
    plt.text(get_x_pos(nh4sh, p), get_y_pos(nh4sh, p), "NH4SH")
    plt.yticks([])
    plt.legend(fontsize=8)
    plt.show()


run_atmosphere_graphs()
