import numpy as np
import matplotlib.pyplot as plt
import unpack as up
import jackknife as jk


class DataContainer:
    def __init__(self):
        self.channels = None

        self.tb = None
        self.tbe_slope = None
        self.tbe_offset = None
        self.frequencies = None
        self.wavelengths = None
        self.dates = None
        self.j_el = None
        self.m_el = None
        self.unpacker_type = up.Unpack

    def get_unpacker(self):
        unpacker = self.unpacker_type().prepare().adjust()
        self.channels = unpacker.channel_obj_list
        return unpacker

    def populate(self):
        unpacker = self.get_unpacker()
        self.tb, self.tbe_slope, self.tbe_offset = unpacker.get_temperatures()
        self.frequencies = unpacker.get_frequencies()
        self.wavelengths = unpacker.get_wavelengths()
        self.dates = unpacker.get_dates()
        self.j_el = unpacker.get_j_el()
        self.m_el = unpacker.get_m_el()
        return self

    def get_original_flux_errors(self):
        return [channel.error for channel in self.channels]

    def yield_model_comp_info(self):
        return self.frequencies, self.tb, self.tbe_slope, self.tbe_offset

    def old_yield_model_comp(self, model_name):
        return self.frequencies, self.tb, self.tbe_slope, self.tbe_offset, model_name
