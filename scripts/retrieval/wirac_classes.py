#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Integration and DataRetrieval classes implementation for WIRAC instrument

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
from xarray.core import dataarray

from gromora_integration import Integration
from gromora_retrievals import DataRetrieval
import GROMORA_library as GROMORA_library

def return_bad_channels_wirac(date):
    '''
    to get the bad channels as a function of the date for WIRAC
    
    Parameters
    ----------
    date : datetime object
        DESCRIPTION.
    
    '''
    #if year == 2019,....
    bad_channels = [8192, 17408, 19456, 23552, 27648, 29696, 31744]
    return bad_channels

class IntegrationWIRAC(Integration):
    '''
    Implementing the Integration class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, 
        date, 
        basename_lvl1, 
        integration_strategy = None,
        integration_time = 1
        ):

        observation_frequency = 1.4217504e11
        instrument_name = "GROMOS"

        self.bandwidth = [1e9]
        spectrometers = ["AC240"]
        
        self.lo = 1.45875e11

        level1_folder = os.path.join(basename_lvl1, instrument_name)

        #integration_strategy = 'simple'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder)
    
    def return_bad_channels(self, date, spectro):
        return return_bad_channels_wirac(date)

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1'):
        '''
        Correction function for the troposphere. 
        
        Invidual correction for each spectrometers specified !
        '''
        return GROMORA_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

class WIRAC_LvL2(DataRetrieval):
    '''
    Implementing the Dataretrieval class for the WIRAC case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, integration_strategy, integration_time, extra_base=''):
        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "WIRAC"

        self.bandwidth = [1e9]
        spectrometers = ["AC240"]

        self.lo = 1.45875e11
        self.reference_elevation_angle = 90
        
        level1_folder = basename_lvl1 # os.path.join(basename_lvl1, instrument_name)
        level2_folder = basename_lvl2# os.path.join(basename_lvl2, instrument_name)

        # Can be used for plotting names (GROMOS_AC240_...)
        self.basename_plot_level2 = instrument_name+'_'+spectrometers[0]+'_'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder, level2_folder, extra_base)
    
    def return_bad_channels(self, date, spectro):
        return return_bad_channels_wirac(date)

    
    def baseline_period(self, retrieval_param):
        '''
        Depending on the dates, function to apply the appropriate baseline periods for the WIRAC retrievals

        ''' 
        # if (retrieval_param['date'] >= datetime.date(2019,1,15)) & (retrieval_param['date'] < datetime.date(2019,1,16)):
        #     baseline_periods = np.array([400e6])
        # else:
        #     baseline_periods = np.array([])
            
        return np.array([])
    
    @property
    def basecolor(self):
       return '#d7191c' 

    @property
    def cost_threshold(self):
        return 0.1

    @property
    def polyfit_threshold(self):
        return 0.1