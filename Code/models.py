import numpy as np

PATH = "/home/ramsey/Documents/Research/Jupiter/SZA/models/"


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
