#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Classes for GROMOS instrument

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
from xarray.core import dataarray

from base_classes import Integration# , DataRetrieval

from gromora_retrievals import DataRetrieval
import GROSOM_library

def return_bad_channels_gromos(date):
    '''
    to get the bad channels as a function of the date for GROMOS
    
    Parameters
    ----------
    date : datetime object
        DESCRIPTION.
    
    '''
    #if year == 2019,....
    bad_channels = np.arange(16383,16384)
    return bad_channels

class IntegrationGROMOS(Integration):
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
        spectrometers = ["AC240"]
        
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
        return GROSOM_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

class GROMOS_LvL2(DataRetrieval):
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
        spectrometers = ["AC240"]

        self.lo = 1.45875e11
        self.reference_elevation_angle = 90
        
        level1_folder = basename_lvl1 # os.path.join(basename_lvl1, instrument_name)
        level2_folder = basename_lvl2# os.path.join(basename_lvl2, instrument_name)

        # Can be used for plotting names (GROMOS_AC240_...)
        self.basename_plot_level2 = instrument_name+'_'+spectrometers[0]+'_'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder, level2_folder, extra_base)
    
    def return_bad_channels(self, date, spectro):

        return return_bad_channels_gromos(date)

    
    def baseline_period(self, retrieval_param):
        '''
        Depending on the dates, function to apply the appropriate baseline periods for the GROMOS retrievals

        ''' 
        if (retrieval_param['date'] >= datetime.date(2019,1,15)) & (retrieval_param['date'] < datetime.date(2019,1,16)):
            baseline_periods = np.array([400e6])

        return baseline_periods

    def make_f_grid_double_sideband(self, retrieval_param): 
        '''
        create simulation frequency grid

        '''
        usb_grid= np.arange(148.975e9,150.175e9,100e6)

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