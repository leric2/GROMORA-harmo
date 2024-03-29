#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Integration and DataRetrieval classes implementation for SOMORA instrument


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

from gromora_integration import Integration
from gromora_retrievals import DataRetrieval
import GROMORA_library as GROMORA_library

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
        return GROMORA_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

class SOMORA_LvL2(DataRetrieval):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, integration_strategy, integration_time, extra_base=''):
        '''
        Some specific parameters to implement for the SOMORA instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "SOMORA"

        self.bandwidth = [1e9]
        spectrometers = ["AC240"]

        self.lo = 1.49275e11
        self.reference_elevation_angle = 90
        self.antenna_fwhm = 2.5
        
        level1_folder = basename_lvl1#  os.path.join(basename_lvl1, instrument_name)
        level2_folder = basename_lvl2#  os.path.join(basename_lvl2, instrument_name)^

        self.institution = 'Federal Office of Meteorology and Climatology;METEOSWISS'
        self.affiliation = 'meteoswiss001'
        self.source = 'MWR.O3_METEOSWISS001_HARMON.2022'
        self.location = 'PAYERNE'
        self.longitude = 6.94
        self.latitude = 46.82
        self.altitude = 491
        self.timezone = 'Z'

        self.name_PI = 'Maillard Barras;Eliane'
        self.contact = 'Ch. de l\'Aerologie;CH-1530 Payerne;SWITZERLAND'

        # Can be used for plotting names (SOMORA_AC240_...)
        self.basename_plot_level2 = instrument_name+'_'+spectrometers[0]+'_'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder, level2_folder, extra_base)

    def return_bad_channels(self, date, spectro):

        return return_bad_channels_somora(date)  
    
    def baseline_period(self, retrieval_param):
        '''
        Depending on the dates, function to apply the appropriate baseline periods for the GROMOS retrievals

        ''' 
        if (retrieval_param['date'] >= datetime.date(2009,1,1)) & (retrieval_param['date'] < datetime.date(2019,10,1)):
            baseline_periods = np.array([110e6, 216e6, 310e6])
        elif (retrieval_param['date'] >= datetime.date(2019,10,1)) & (retrieval_param['date'] < datetime.date(2020,4,15)):
            baseline_periods = np.array([110e6, 186e6, 347e6, 888e6]) #888e6
        elif (retrieval_param['date'] >= datetime.date(2020,4,15)) & (retrieval_param['date'] < datetime.date(2020,9,30)):
            baseline_periods = np.array([92e6, 110e6, 364e6, 400e6]) #baseline_periods = np.array([110e6, 186e6, 364e6])
        elif (retrieval_param['date'] >= datetime.date(2020,9,30)):
            baseline_periods = np.array([96e6, 110e6, 186e6, 347e6, 382e6]) # 96e6 888e6 baseline_periods = np.array([110e6, 186e6, 364e6])
        else:
            baseline_periods = np.array([])

        return baseline_periods

    def correct_pointing(self, retrieval_param):
        return 0

    def make_f_grid_double_sideband(self, retrieval_param): 
        '''
        Create simulation frequency grid when the sideband response is included.

        '''
        usb_grid = self.usb_grid
        n_f = retrieval_param["number_of_freq_points"]  # Number of points
        bw = 1.3*retrieval_param["bandwidth"]  # Bandwidth
        x = np.linspace(-1, 1, n_f)
        f_grid = x ** 3 + x / retrieval_param["irregularity_f_grid"]
        f_grid = f_grid * bw / (max(f_grid) - min(f_grid)) + \
            retrieval_param['obs_freq']

        #f_grid = np.linspace(retrieval_param["f_min"]-10, retrieval_param["f_max"]+10, n_f)
        f_grid = np.concatenate((f_grid, usb_grid))
        if retrieval_param["show_f_grid"]:
            fig = plt.figure()
            plt.semilogy(f_grid[1:]/1e9, np.diff(f_grid)/1e3, '.')
            # plt.xlim((retrieval_param['obs_freq']-200e6) /
            #          1e9, (retrieval_param['obs_freq']+200e6)/1e9)
            # plt.ylim(0,300)
            plt.ylabel(r'$\Delta f$ [kHz]')
            plt.suptitle('Frequency grid spacing')
            plt.show()
        return f_grid

    def cost_threshold(self, year):
        '''
        Cost threshold over which we flag the level 2
        '''
        if year > 2010:
            return 0.1 #0.5 for 2010
        else:
            return 0.5

    @property
    def day2flag_level2(self):
        '''
        A selection of days to flags for the level2 GROMOS data. 
        These days have been identified in the GROMORA time series detailed analysis that can be found in the GROMORA retrievals UG.

        '''
        date2flag =  [
             slice('2012-04-24','2012-04-27')
        ]
        date2flag.append(slice('2016-09-29','2016-11-03'))
        date2flag.append(slice('2018-01-31','2018-02-11'))
        date2flag.append(slice('2018-07-25','2018-08-23'))
        return date2flag

    @property
    def basecolor(self):
       return '#2c7bb6'

    @property
    def polyfit_threshold(self):
        return 0.1

    @property
    def usb_grid(self):
        return np.arange(155.875e9,157.075e9,100e6)

    @property
    def standard_air_pressure(self):
        return 960

    @property
    def standard_air_temperature(self):
        return 10

    @property
    def cycle_duration(self):
        return 2/3600

    @property
    def global_attributes_ndacc(self):
        pi_name='Maillard Barras;Eliane'
        pi_mail='eliane.maillard@meteoswiss.ch'
        do_name = 'Haefele;Alexander'
        do_mail = 'alexander.haefele@meteoswiss.ch'
        rou= 'Please contact Eliane Maillard Barras at eliane.maillard@meteoswiss.ch'
        ackn ='The ozone microwave radiometer SOMORA is operated by MeteoSwiss, Switzerland.'
        description='Atmospheric ozone profiles from continuous measurements by ground-based 142 GHz microwave radiometer SOMORA, Payerne, Switzerland'
        contact = "Ch. de l\'Aerologie;CH-1530 Payerne;SWITZERLAND"
        return dict(
            PI_NAME=pi_name,
            PI_AFFILIATION=self.institution,
            PI_ADDRESS = contact,
            PI_EMAIL = pi_mail,
            DO_NAME= do_name,
            DO_AFFILIATION=self.institution,
            DO_ADDRESS = contact,
            DO_EMAIL =  do_mail,
            DS_NAME= pi_name,
            DS_AFFILIATION=self.institution,
            DS_ADDRESS = contact,
            DS_EMAIL = pi_mail,
            DATA_DESCRIPTION = description,
            DATA_RULES_OF_USE =rou,
            DATA_ACKNOWLEDGEMENT=ackn,
            FILE_PROJECT_ID = 'NDACC-SOMORA'
        )
    # def define_retrieval_param(self, retrieval_param):
    #     '''
    #     Parameters
    #     ----------
    #     retrieval_param : dict
    #         This function fills the dict containing all parameters for GROMOS retrievals.

    #     Returns
    #     -------
    #     retrieval_param

    #     '''
    #     return self.define_retrieval_param

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
  