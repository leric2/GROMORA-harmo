#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Classes for SOMORA instrument

Example:
    E...

        $ python example_google.py

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo: all

"""
from abc import ABC
import os
import datetime

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from base_classes import Integration, DataRetrieval
import GROSOM_library

def return_bad_channels_somora(date):
    '''
    to get the bad channels as a function of the date for SOMORA
    
    Parameters
    ----------
    date : datetime object
        DESCRIPTION.
    
    '''
    #if year == 2019,....
    bad_channels = np.arange(0,104)
    return bad_channels

class IntegrationSOMORA(Integration):
    '''
    Implementing the Integration class for the SOMORA case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, 
        date, 
        basename_lvl1, 
        integration_strategy = None,
        integration_time = 1
        ):

        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "SOMORA"

        self.bandwidth = [1e9]
        spectrometers = ["AC240"]
        
        level1_folder = os.path.join(basename_lvl1, instrument_name)

        #integration_strategy = 'simple'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder)
    
    def return_bad_channels(self, date, spectro):

        return return_bad_channels_somora(date)

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1'):
        '''
        Correction function for the troposphere. 
        
        Invidual correction for each spectrometers specified !
        '''
        return GROSOM_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

class SOMORA_LvL2(DataRetrieval):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, integration_strategy, integration_time):
        '''
        Some specific parameters to implement for the SOMORA instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "SOMORA"

        self.bandwidth = [1e9]
        spectrometers = ["AC240"]

        self.lo = 1.49275e11
        
        level1_folder = basename_lvl1#  os.path.join(basename_lvl1, instrument_name)
        level2_folder = basename_lvl2#  os.path.join(basename_lvl2, instrument_name)

        # Can be used for plotting names (SOMORA_AC240_...)
        self.basename_plot_level2 = instrument_name+'_'+spectrometers[0]+'_'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder, level2_folder)

    def return_bad_channels(self, date, spectro):

        return return_bad_channels_somora(date)  
    
    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1'):
        '''
        Correction function for the troposphere. 
        
        Invidual correction for each spectrometers specified !
        '''
        return GROSOM_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

    # def find_bad_channels(self, Tb_min, Tb_max, boxcar_size, boxcar_thresh):
    #     '''
    #     Parameters
    #     ----------
    #     level1b_dataset : TYPE
    #         DESCRIPTION.

    #     Returns
    #     -------
    #     None.

    #     '''
    #     bad_channels = self.return_bad_channel_SOMORA(self.date)
    #     #self.level1b_ds = super().find_bad_channels(bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
    #     #
    #     return data_GROSOM.find_bad_channels(self.level1b_ds, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
  