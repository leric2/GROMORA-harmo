#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Integration and DataRetrieval classes implementation for GROMOS instrument


"""
from abc import ABC
import os, datetime
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

def return_bad_channels_gromos(date):
    '''
    to get the bad channels as a function of the date for GROMOS
    
    Parameters
    ----------
    date : datetime object
    
    '''
    
    bad_channels = np.arange(16383,16384)
    return bad_channels

class IntegrationGROMOSFB(Integration):
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

        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "GROMOS"

        self.bandwidth = [1e9]
        spectrometers = ["FB"]
        
        self.lo = 1.45875e11

        level1_folder = os.path.join(basename_lvl1, instrument_name)

        #integration_strategy = 'simple'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder)
    
    def return_bad_channels(self, date, spectro):
        return return_bad_channels_gromos(date)

    # def compare_Tb_chunks(self, dim='time', idx=[0], save = False, Tb_corr = False):
    #     figures = list()

    #     for i in idx:
    #         figures.append(GROSOM_library.plot_Tb_chunks(self, self.integrated_dataset, i)) 

    #     if T_corr:
    #         for i in idx:
    #             figures.append(GROSOM_library.plot_Tb_corr_chunks(self, self.integrated_dataset, i)) 
        
    #     if save:
    #         raise NotImplementedError()

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1'):
        '''
        Correction function for the troposphere. 
        
        Invidual correction for each spectrometers specified !
        '''
        return GROMORA_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

class GROMOS_FB_LvL2(DataRetrieval):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, integration_strategy, integration_time, extra_base=''):
        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        observation_frequency = 1.4217504e11
        instrument_name = "GROMOS"

        self.bandwidth = [1e9]
        spectrometers = ["FB"]

        
        self.lo = 1.45875e11
        self.reference_elevation_angle = 90
        
        level1_folder = basename_lvl1 # os.path.join(basename_lvl1, instrument_name)
        level2_folder = basename_lvl2# os.path.join(basename_lvl2, instrument_name)
        
        self.institution = 'University of Bern;UBERN'
        self.affiliation = 'ubern001'
        self.source = 'MWR.O3_UBERN'
        self.longitude = 7.44
        self.latitude = 46.95
        self.altitude = 560
        
        # Can be used for plotting names (GROMOS_AC240_...)
        self.basename_plot_level2 = instrument_name+'_'+spectrometers[0]+'_'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder, level2_folder, extra_base)
    
    def return_bad_channels(self, date, spectro):
        return return_bad_channels_gromos(date)

    
    def baseline_period(self, retrieval_param):
        '''
        Depending on the dates, function to apply the appropriate baseline periods for the GROMOS retrievals

        ''' 
        #if (retrieval_param['date'] >= datetime.date(2009,1,1)) & (retrieval_param['date'] < datetime.date(2015,2,23)):
        #     baseline_periods = np.array([178e6, 240e6, 360e6])
        # elif  (retrieval_param['date'] >= datetime.date(2015,2,23)) & (retrieval_param['date'] < datetime.date(2015,8,31)):
        #     baseline_periods = np.array([140e6, 240e6, 400e6])
        # elif  (retrieval_param['date'] >= datetime.date(2015,8,31)) & (retrieval_param['date'] < datetime.date(2017,1,1)):
        #     baseline_periods = np.array([160e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2017,1,1)) & (retrieval_param['date'] < datetime.date(2018,1,1)):
        #     baseline_periods = np.array([178e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2018,1,1)) & (retrieval_param['date'] < datetime.date(2019,1,1)):
        #     baseline_periods = np.array([135e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2019,1,1)) & (retrieval_param['date'] < datetime.date(2019,3,15)):
        #     baseline_periods = np.array([155e6, 240e6, 360e6])
        # elif (retrieval_param['date'] >= datetime.date(2019,3,15)) & (retrieval_param['date'] < datetime.date(2022,2,16)):
        #     baseline_periods = np.array([135e6, 178e6, 240e6])
        # elif retrieval_param['date'] > datetime.date(2022,2,16):
        #     baseline_periods = np.array([135e6, 178e6, 240e6])
        #else:
        baseline_periods = np.array([])
            #raise ValueError('No Sinefit implement for FB yet')
            
        return baseline_periods
    
    def correct_pointing(self, retrieval_param):
        if (retrieval_param['date'] >= datetime.date(2019,2,12)) & (retrieval_param['date'] < datetime.date(2019,3,13)):
            return -5
        else:
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
        return 0.1 

    # @property
    # def day2flag_level2(self):
    #     '''
    #     A selection of days to flags for the level2 GROMOS data. 
    #     These days have been identified in the GROMORA time series detailed analysis that can be found in the GROMORA retrievals UG.

    #     '''
    #     date2flag_gromos =  [
    #         datetime.date(2015,8,26), datetime.date(2015,8,27), datetime.date(2015,8,28),
    #         pd.date_range('2012-07-24', '2012-08-07'),
    #         pd.date_range('2019-01-14', '2019-02-12'),

    #     ]
    #     return date2flag_gromos

    @property
    def usb_grid(self):
        return np.arange(148.975e9,150.175e9,100e6)

    @property
    def basecolor(self):
       return '#d7191c' 

    @property
    def polyfit_threshold(self):
        return 0.1

    @property
    def standard_air_pressure(self):
        return 955
        
    @property
    def standard_air_temperature(self):
        return 10

    @property
    def cycle_duration(self):
        return 7/3600
    
    @property
    def global_attributes_ndacc(self):
        pi_mail='axel.murk@unibe.ch'
        do_name = 'Sauvageat;Eric'
        do_mail = 'eric.sauvageat@unibe.ch'
        rou= 'Please contact Axel Murk at axel.murk@unibe.ch'
        ackn ='The middle atmospheric ozone radiometer GROMOS is operated by the Institute of Applied Physics, University of Bern, Switzerland.'
        description='Atmospheric ozone profiles from continuous measurements of ground-based 142 GHz-microwave radiometer GROMOS at Bern, Switzerland'
        return dict(
            PI_NAME=self.name_PI,
            PI_AFFILIATION=self.institution,
            PI_ADDRESS = self.contact,
            PI_EMAIL = pi_mail,
            DO_NAME= do_name,
            DO_AFFILIATION=self.institution,
            DO_ADDRESS = self.contact,
            DO_EMAIL =  do_mail,
            DS_NAME= do_name,
            DS_AFFILIATION=self.institution,
            DS_ADDRESS = self.contact,
            DS_EMAIL = do_mail,
            DATA_DESCRIPTION = description,
            DATA_RULES_OF_USE =rou,
            DATA_ACKNOWLEDGEMENT=ackn,
            FILE_PROJECT_ID = 'NDACC-GROMOS'
        )