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
from utils_GROSOM import save_pngs
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
            #figures.append(mopi5_library.compare_Tb_mopi5(self, self.integrated_data, i)) 
            if with_corr:
                figures.append(mopi5_library.compare_spectra_mopi5_new(self, self.integrated_data, i, title=title))
            else:
                figures.append(mopi5_library.compare_spectra_only_mopi5(self, self.integrated_data, i, title=title))

        if save_plot:
            save_single_pdf(self.level1_folder+'spectra_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)

    def compare_interpolated_spectra_mopi5(self, dim='time', idx=[0], spectrometers=['AC240','USRP-A'], use_basis='U5303', save_plot = False, identifier=[]):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        for i in idx:
            title = self.integration_strategy + ' n.'+ str(i) + ' : ' + str(identifier[i])
            #figures.append(mopi5_library.compare_Tb_mopi5(self, self.integrated_data, i)) 
            figures.append(mopi5_library.compare_spectra_diff_mopi5(self, self.integrated_data, i, spectrometers, use_basis, title=title))

        if save_plot:
            save_single_pdf(self.level1_folder+'spectra_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1', basis_spectro='AC240', skip_ch=[1000,1000], num_of_ch=500, interp=False):
        '''
        Correction function for the troposphere. 
        
        For MOPI, we call it a first time 
        '''
        if interp:
            return mopi5_library.correct_troposphere_interpolated(self, spectrometers, dim, method='Ingold_v1', use_basis=basis_spectro, skip_ch=skip_ch, num_of_ch=num_of_ch)
        else:
            return mopi5_library.correct_troposphere(self, spectrometers, dim, method='Ingold_v1', use_basis=basis_spectro, skip_ch=skip_ch, num_of_ch=num_of_ch)

    def plot_time_min_comp(self):
        return mopi5_library.plot_time_min_comp(self)

    def plot_time_series_all_mopi5(self, title='', special=False):
        figures = list()
        if special:
            figures.append(mopi5_library.plot_ts_mopi5_Feb(self, title=title))
        else:
            figures.append(mopi5_library.plot_ts_mopi5(self, title=title))

        save_single_pdf(self.level1_folder+'time_series_'+self.datestr+'.pdf', figures)
        save_pngs(self.level1_folder+'time_series_'+self.datestr+'_', figures)


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

    def import_standard_retrieval_params_mopi5(self):

        retrieval_param = dict()
        retrieval_param["retrieval_type"] = 1

        retrieval_param["obs_freq"] = self.observation_frequency
    
        retrieval_param["plot_meteo_ds"] = True
        retrieval_param["number_of_freq_points"] = 601
        retrieval_param["irregularity_f_grid"] = 10
        retrieval_param["show_f_grid"] = True

        retrieval_param["z_top_sim_grid"] = 97e3
        retrieval_param["z_bottom_sim_grid"] = 1000
        retrieval_param["z_resolution_sim_grid"] = 1e3

        retrieval_param["z_top_ret_grid"] = 95e3
        retrieval_param["z_bottom_ret_grid"] = 1000
        retrieval_param["z_resolution_ret_grid"] = 3e3

        #retrieval_param["z_top_ret_grid_h2o"] = 50e3
        #retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

        
        retrieval_param['unit_var_y']  = 3


        retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
        retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
        retrieval_param["apriori_O3_cov"] = 1.5e-6

        retrieval_param['water_vapor_model'] = "H2O-PWR98"
        #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
        #retrieval_param["azimuth_angle"]=32

        
        
        retrieval_param['ecmwf_store_location'] ='/scratch/ECMWF'
        retrieval_param['extra_time_ecmwf'] = 6


        return retrieval_param

    def return_bad_channels(self, date, spectro):
        '''
        to get the bad channels as a function of the date for MOPI5
        
        Parameters
        ----------
        date : datetime object
            DESCRIPTION.
        
        '''
        number_of_channel = len(self.integrated_data[spectro].channel_idx)
        #intermediate_frequency = self.calibrated_data[spectro].intermediate_frequency.data
        return mopi5_library.return_bad_channels_mopi5(number_of_channel, date, spectro)

    def compare_spectra_mopi5(self, spectrometers=[], dim='time', idx=[0], save_plot = False, identifier=[], lowerBound=[], with_corr = True, corr_band=[], title=None):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        for i in idx:
            if title is None:
                title = self.integration_strategy + ' n.'+ str(i) + ' : ' + str(identifier[i])+', int_time='+str(10*self.integrated_data['U5303'].chunk_size[i].data) + 'min'
            #figures.append(mopi5_library.compare_Tb_mopi5(self, self.integrated_data, i)) 
            if with_corr:
                figures.append(mopi5_library.compare_spectra_mopi5_new(self, self.integrated_data, spectrometers, i, title=title, corr_band=corr_band))
            else:
                #print(corr_band)
                title ='Integrated spectra with $T_{B,mean}$ between '+str(lowerBound[i])+ ' and '+str(identifier[i])+'K'
                figures.append(mopi5_library.compare_spectra_only_mopi5(self, self.integrated_data, spectrometers, i, title=title, corr_band=corr_band))

        if save_plot:
            save_single_pdf(self.level1_folder+'spectra_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)

    def compare_binned_spectra_mopi5(self, dim='time', idx=[0], save_plot = False, identifier=[], use_basis = 'U5303', df=200e3):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        for i in idx:
            title = self.integration_strategy + ' n.'+ str(i) + ' : ' + str(identifier[i]) + ', binned with df='+str(df/1e3)+' kHz'
            figures.append(mopi5_library.compare_binned_spectra_only_mopi5(self, self.integrated_data, i, use_basis=use_basis, title=title))

        if save_plot:
            save_single_pdf(self.level1_folder+'spectra_binned_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)
    
    def compare_spectra_binned_interp_mopi5_factor(self, dim='time', idx=[0], save_plot = False, spectrometers=['AC240'], identifier=[], use_basis = 'U5303', alpha=[0,7,8,9], binning=8, lowerBound=[], corr_band=[]):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        
        for i in idx:            
            title ='Integrated spectra with $T_{B,mean}$ between '+str(lowerBound[i])+ ' and '+str(identifier[i])+'K'
            figures.append(mopi5_library.compare_spectra_binned_interp_mopi5_clean_factor(self, self.integrated_data, i, spectrometers=spectrometers, use_basis=use_basis, alpha=alpha, binning=binning, title=title, corr_band=corr_band))

        if save_plot:
            save_single_pdf(self.level1_folder+'binned_factor'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)
    
    def compare_spectra_binned_interp_mopi5(self, dim='time', idx=[0], spectrometers=['AC240','USRP-A'], use_basis='U5303', save_plot = False, identifier=[], corrected=False, clean=False):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        for i in idx:
            title = self.integration_strategy + ' n.'+ str(i) + ' : ' + str(identifier[i]) + ', binned and interpolated Tb with ' + use_basis
            #figures.append(mopi5_library.compare_Tb_mopi5(self, self.integrated_data, i)) 
            if clean:
                if corrected:
                    lowerBound = [0, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 190]
                    title ='Integrated and corrected spectra with $T_{B,mean}$ between '+str(lowerBound[i])+ ' and '+str(identifier[i])+'K'
                    figures.append(mopi5_library.compare_spectra_binned_interp_mopi5_clean_corr(self, self.integrated_data, i, spectrometers, use_basis, title=title))
                else:
                    lowerBound = [0, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 190]
                    title ='Integrated spectra with $T_{B,mean}$ between '+str(lowerBound[i])+ ' and '+str(identifier[i])+'K'
                    figures.append(mopi5_library.compare_spectra_binned_interp_mopi5_clean(self, self.integrated_data, i, spectrometers, use_basis, title=title))
            else:
                if corrected:
                    figures.append(mopi5_library.compare_spectra_binned_interp_mopi5_corrected(self, self.integrated_data, i, spectrometers, use_basis, title=title))
                else:
                    figures.append(mopi5_library.compare_spectra_binned_interp_mopi5(self, self.integrated_data, i, spectrometers, use_basis, title=title))
            
        if save_plot:
            if clean:
                if corrected:
                    print('no png saved')
                #else:
                    #save_pngs(self.level1_folder+'spectra_interp_diff_comparison_corr_'+self.integration_strategy+'_'+self.datestr+'_', figures)
            if corrected:
                save_single_pdf(self.level1_folder+'spectra_interp_diff_comparison_corr_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)
            else: 
                save_single_pdf(self.level1_folder+'spectra_interp_diff_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)
    def compare_interpolated_spectra_mopi5(self, dim='time', idx=[0], spectrometers=['AC240','USRP-A'], use_basis='U5303', save_plot = False, identifier=[]):
        figures = list()
        #spectro_ds = self.calibrated_data[s]
        for i in idx:
            title = self.integration_strategy + ' n.'+ str(i) + ' : ' + str(identifier[i]) + ', interpolated Tb with ' + use_basis
            #figures.append(mopi5_library.compare_Tb_mopi5(self, self.integrated_data, i)) 
            figures.append(mopi5_library.compare_spectra_diff_mopi5(self, self.integrated_data, i, spectrometers, use_basis, title=title))

        if save_plot:
            save_single_pdf(self.level1_folder+'spectra_interp_diff_comparison_'+self.integration_strategy+'_'+self.datestr+'_'+str(idx)+'.pdf', figures)


    def correction_function_mopi5(self, spectro_as_basis='U5303', t_trop=290):
        '''
        From Jonas
        Generate a correction function, given frequencies, brightness temperatures
        and a weighted mean torpospheric temperature.
        Returns a function f: (f, y) -> y_corr
        '''
        l = len(self.integrated_data[spectro_as_basis].time)
        y_corr = dict()
        for i, s in enumerate(["U5303", "AC240", "USRP-A"]):
            y_corr[s] = np.ones((l,len(self.integrated_data[s].channel_idx.data)))*np.nan

        for c in range(l):
            good_channels = self.integrated_data[spectro_as_basis].good_channels[c].data == 1
            f_basis = self.integrated_data[spectro_as_basis].frequencies.data[good_channels]
            y_basis = self.integrated_data[spectro_as_basis].Tb[c].data[good_channels]

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

                f = self.integrated_data[s].frequencies.data
                y = self.integrated_data[s].Tb[c].data

                # Correct by frequency dependant opacity for the three spectrometers
                y_eff = a * f + b

                exp_neg_tau = (t_trop - y_eff) / (t_trop - 2.7)
                # tau = -np.log(exp_neg_tau)
                y_corr[s][c,:]  = (y - t_trop * (1 - exp_neg_tau)) / exp_neg_tau

        for i, spectro in enumerate(["U5303", "AC240", "USRP-A"]):
            new_tb_corr = xr.DataArray(y_corr[spectro], dims = ['time', 'channel_idx'])
            self.integrated_data[spectro] = self.integrated_data[spectro].assign(Tb_corr_old = self.integrated_data[spectro].Tb_corr.rename('old_tb_corr'))
            self.integrated_data[spectro] = self.integrated_data[spectro].assign(Tb_corr = new_tb_corr)

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
        return mopi5_retrievals.retrieve_cycle_tropospheric_corrected_mopi5(spectro_dataset, retrieval_param)

    def plot_level2_from_tropospheric_corrected_spectra(self, ac, spectro_dataset, retrieval_param, title, figure_list, fshift, bl):
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        ac : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.
        title : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        return mopi5_library.plot_level2_from_tropospheric_corrected_mopi5(spectro_dataset, ac, retrieval_param, title, figure_list, fshift=fshift, compute_bl=bl)

    def plot_time_min_comp(self):
        return mopi5_library.plot_time_min_comp(self)
