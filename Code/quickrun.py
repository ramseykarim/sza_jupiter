"""
Quick run script for plotting T_B data
"""
import unpack as up
import plotting as pltt
import matplotlib.pyplot as plt
import error_investigation as ei


def run_error_comp():
    p = ei.FlexPlotting()
    u_bootstrap = up.Unpack().error_investigation_init()
    results = u_bootstrap.error_investigation.execute_flux_testing(100)
    p.add_flux_scatter(results)
    unpackers = {}
    for i in range(4):
        u = ei.FlexUnpack().assign(i+1).adjust()
        unpackers[i+1] = u
        p.add_my_data(u)
    p.show_prep()
    ei.add_best_mean_bars(unpackers[1], results)
    p.show()


def run_normal():
    u = up.Unpack().adjust()
    p = pltt.Plotting()
    p.add_my_data(u)
    p.add_all_other_data()
    p.add_model_plot("saturated", color='k')
    p.add_model_plot("50abv2bar", color='green')
    u.print_points()
    p.show()


def run_model_comp():
    p = pltt.Plotting()
    plt.figure(1)
    plt.subplot(121)
    p.add_model_plot("saturated", color='k')
    p.add_model_plot("50abv1bar", color='blue')
    p.add_model_plot("50abv2bar", color='red')
    p.add_model_plot("80abv2bar", color='green')
    p.add_model_plot("saturated_newnh3", color='k', line_sty='-.')
    p.add_model_plot("50abv1bar_newnh3", color='blue', line_sty='-.')
    p.add_model_plot("50abv2bar_newnh3", color='red', line_sty='-.')
    p.add_model_plot("80abv2bar_newnh3", color='green', line_sty='-.')
    plt.title("Dotted - new nh3, Solid - old nh3")
    plt.subplot(122)
    p.add_model_compare("saturated", color='k')
    p.add_model_compare("50abv1bar", color='blue')
    p.add_model_compare("50abv2bar", color='red')
    p.add_model_compare("80abv2bar", color='green')
    plt.title("New nh3 - old nh3")
    p.show_prep()
    p.show()

run_normal()
