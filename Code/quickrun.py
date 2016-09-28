"""
Quick run script for plotting T_B data
"""
import unpack as up
import plotting as pltt

u = up.Unpack()
p = pltt.Plotting()
p.add_my_data(u)
p.add_all_other_data()
p.add_model_plot("saturated", color='k')
p.show()
