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

from base_classes import Integration, DataRetrieval
import mopi5_retrievals
import mopi5_library
from utils_GROSOM import save_single_pdf
import retrieval_module

class IntegrationMOPI5(Integration):
    '''
    Implementing the Integration class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, 
        date, 
        basename_lvl1, 
        integration_strategy = None,
        integration_time = 1,
        spectrometers = []
        ):

        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        observation_frequency = 110.836e9
        instrument_name = "mopi5"

        self.bandwidth = [1.6e9, 1e9, 200e6]
        #self.number_of_channel = [16384, 16384, 16384]
        if not spectrometers:
            spectrometers = ["U5303","AC240","USRP-A"]
        
        level1_folder = basename_lvl1

        #integration_strategy = 'simple'

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder)


    def return_bad_channels(self, date, spectro):
        '''
        to get the bad channels as a function of the date for MOPI5
        
        Parameters
        ----------
        date : datetime object
            DESCRIPTION.
        
        '''
        number_of_channel = len(self.calibrated_data[spectro].channel_idx)
        #intermediate_frequency = self.calibrated_data[spectro].intermediate_frequency.data
        return mopi5_library.return_bad_channels_mopi5(number_of_channel, date, spectro)

    def compare_spectra_mopi5(self, dim='time', idx=[0], save_plot = False, identifier=[], with_corr = True):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        for i in idx:
            title = self.integration_strategy + ' n.'+ str(i) + ' : ' + str(identifier[i])
            #figures.append(mopi5_library.compare_Tb_mopi5(self, self.integrated_dataset, i)) 
            if with_corr:
                figures.append(mopi5_library.compare_spectra_mopi5_new(self, self.integrated_dataset, i, title=title))
            else:
                figures.append(mopi5_library.compare_spectra_only_mopi5(self, self.integrated_dataset, i, title=title))

        if save_plot:
            save_single_pdf(self.level1_folder+'spectra_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)
    
    def correct_troposphere_old(self, spectrometers, dim, method='Ingold_v1'):
        '''
        Correction function for the troposphere. 
        
        Invidual correction for each spectrometers specified !
        '''
        raise NotImplementedError
        import GROSOM_library
        print('Using the correction function from GROSOM, TODO')
        return GROSOM_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1')

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1', basis_spectro='AC240'):
        '''
        Correction function for the troposphere. 
        
        For MOPI, we call it a first time 
        '''
        #basis_spectro = 'AC240'
        return mopi5_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1', use_basis=basis_spectro)

class MOPI5_LvL2(DataRetrieval):
    '''
    Implementing the Dataretrieval class for the MOPI5 case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, integration_strategy, integration_time):
        '''
        Some specific parameters to implement for the SOMORA instances (only constant stuff...)
        '''
        observation_frequency = 110.836e9
        instrument_name = "mopi5"

        self.bandwidth = [1.6e9, 1e9, 200e6]
        #self.number_of_channel = [16384, 16384, 16384]
        spectrometers = ["U5303","AC240","USRP-A"]
        
        level1_folder = basename_lvl1
        level2_folder = basename_lvl2

        super().__init__(instrument_name, observation_frequency, spectrometers, integration_strategy, integration_time, date, level1_folder, level2_folder)

    def return_bad_channels(self, date, spectro):
        '''
        to get the bad channels as a function of the date for MOPI5
        
        Parameters
        ----------
        date : datetime object
            DESCRIPTION.
        
        '''
        number_of_channel = len(self.integrated_dataset[spectro].channel_idx)
        #intermediate_frequency = self.calibrated_data[spectro].intermediate_frequency.data
        return mopi5_library.return_bad_channels_mopi5(number_of_channel, date, spectro)

    def compare_spectra_mopi5(self, dim='time', idx=[0], save_plot = False, identifier=[], with_corr = True):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        for i in idx:
            title = self.integration_strategy + ' n.'+ str(i) + ' : ' + str(identifier[i])
            #figures.append(mopi5_library.compare_Tb_mopi5(self, self.integrated_dataset, i)) 
            if with_corr:
                figures.append(mopi5_library.compare_spectra_mopi5_new(self, self.integrated_dataset, i, title=title))
            else:
                figures.append(mopi5_library.compare_spectra_only_mopi5(self, self.integrated_dataset, i, title=title))

        if save_plot:
            save_single_pdf(self.level1_folder+'spectra_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)
    
    def correction_function_mopi5(self, spectro_as_basis='U5303', t_trop=290):
        '''
        From Jonas
        Generate a correction function, given frequencies, brightness temperatures
        and a weighted mean torpospheric temperature.
        Returns a function f: (f, y) -> y_corr
        '''
        l = len(self.data[spectro_as_basis].time)
        y_corr = dict()
        for i, s in enumerate(["U5303", "AC240", "USRP-A"]):
            y_corr[s] = np.ones((l,len(self.data[s].channel_idx.data)))*np.nan

        for c in range(l):
            good_channels = self.data[spectro_as_basis].good_channels[c].data == 1
            f_basis = self.data[spectro_as_basis].frequencies.data[good_channels]
            y_basis = self.data[spectro_as_basis].Tb[c].data[good_channels]

            # Skip 100 channels at beginning
            # use 500 channels for reference
            wings = (slice(100, 600), slice(-600, -100))

            yw = np.stack([y_basis[s].mean() for s in wings])
            fw = np.stack([f_basis[s].mean() for s in wings])

            # Fit a line with slope to the spectrum
            # y = a*y0 + b
            a = (yw[1] - yw[0]) / (fw[1] - fw[0])
            b = yw[1] - a * fw[1]

            for i, s in enumerate(["U5303", "AC240", "USRP-A"]):

                f = self.data[s].frequencies.data
                y = self.data[s].Tb[c].data

                # Correct by frequency dependant opacity for the three spectrometers
                y_eff = a * f + b

                exp_neg_tau = (t_trop - y_eff) / (t_trop - 2.7)
                # tau = -np.log(exp_neg_tau)
                y_corr[s][c,:]  = (y - t_trop * (1 - exp_neg_tau)) / exp_neg_tau

        for i, spectro in enumerate(["U5303", "AC240", "USRP-A"]):
            new_tb_corr = xr.DataArray(y_corr[spectro], dims = ['time', 'channel_idx'])
            self.data[spectro] = self.data[spectro].assign(Tb_corr_old = self.data[spectro].Tb_corr.rename('old_tb_corr'))
            self.data[spectro] = self.data[spectro].assign(Tb_corr = new_tb_corr)

        return self

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1', basis_spectro='AC240'):
        '''
        Correction function for the troposphere. 
        
        For MOPI, we call it a first time 
        '''
        #basis_spectro = 'AC240'
        return mopi5_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1', use_basis=basis_spectro)

    def retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''
        retrieval_param["binned_ch"] = False
        retrieval_param['ref_elevation_angle'] = 180
        return retrieval_module.retrieve_cycle_tropospheric_corrected_mopi5(spectro_dataset, retrieval_param)