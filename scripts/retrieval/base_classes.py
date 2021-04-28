#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Main classes definition for radiometer retrieval

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
from utils_GROSOM import save_single_pdf

from dotenv import load_dotenv

# For ARTS, we need to specify some paths
#load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

import GROSOM_library
import retrieval_module
import mopi5_retrievals

from retrievals import data

#%%

# class MicrowaveRadiometerIAP(ABC):
#     '''
#     The base abstract class for every IAP MW radiometer. 

#     It has multiple children classes, for integration, retrieval, etc...

#     '''
#     def __init__(
#         self,
#         instrument_name=None,
#         observation_frequency=None,
#         spectrometers=None,
#         ):
        
#         self.instrument_name = instrument_name
#         self.observation_frequency = observation_frequency
#         self.spectrometers = spectrometers

class Integration(ABC):
    '''
    create another class for integration of level1a data
    
    '''
    def __init__(
        self,
        instrument_name=None,
        observation_frequency=None,
        spectrometers=None,
        integration_strategy='classic',
        integration_time=1,
        dates=None,
        level1_folder=None,
        ):

        self.instrument_name = instrument_name
        self.observation_frequency = observation_frequency
        self.spectrometers = spectrometers
        self.integration_strategy = integration_strategy
        self.int_time = integration_time
        self.level1_folder = level1_folder

        #super().__init__(instrument_name, observation_frequency, spectrometers)

        self.dates = dates

        try:
            len(self.dates) > 1
        except TypeError:
            self.multiday = False
        else:
            self.multiday = True

        # if len(self.dates) > 1:
        #     self.multiday = True
        # else:
        #     self.multiday = False

        self.calibrated_data = dict()
        self.meteo_data = dict()
        self.calibration_flags = dict()
        self.filename_level1a = dict()
        self.filename_level1b = dict()

        for s in self.spectrometers:
            self.filename_level1a[s] = list()
            if self.multiday:
                for date in dates:
                    self.datestr = date.strftime('%Y_%m_%d')
                    self.filename_level1a[s].append(
                        os.path.join(
                    self.level1_folder,
                    self.instrument_name + "_level1a_" +
                    s + "_" + self.datestr
                    ))

                    if self.integration_strategy == 'classic':
                        if self.int_time == 1:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_" +
                            s + "_" + self.datestr
                            )
                        else:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                            s + "_" + self.datestr
                            )
                    else:
                        self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                            s + "_" + self.datestr
                            )
            else:
                self.datestr = self.dates.strftime('%Y_%m_%d')
                self.filename_level1a[s].append(
                    os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1a_" +
                s + "_" + self.datestr
                ))
                if self.integration_strategy == 'classic':
                    if self.int_time == 1:
                        self.filename_level1b[s] = os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_" +
                        s + "_" + self.datestr
                        )
                    else:
                        self.filename_level1b[s] = os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                        s + "_" + self.datestr
                        )
                else:
                    self.filename_level1b[s] = os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                        s + "_" + self.datestr
                        )

    def read_level1a(self):
        ''' 
        Reading level1a in netCDF format
        
        '''
            
        for s in self.spectrometers:
            try:
                #self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = GROSOM_library.read_level1(self._filename_lvl1b)
                self.calibrated_data[s], self.calibration_flags[s], self.meteo_data[s], global_attrs_level1a = GROSOM_library.read_level1(self.filename_level1a[s][0])
                print('Read : ', self.filename_level1a[s][0])
            except FileNotFoundError:
                print('No file for this day, skipping ', self.filename_level1a[s][0])
                
        # Meta data
        self.institution = global_attrs_level1a['institution']
        self.instrument_name_from_level1a = global_attrs_level1a['instrument']
        self.location = global_attrs_level1a['location']

        self.number_of_spectrometer = global_attrs_level1a['number_of_spectrometer']
        self.calibration_version = global_attrs_level1a['calibration_version']

        self.raw_data_filename = global_attrs_level1a['raw_filename']
        self.raw_data_software_version = global_attrs_level1a['raw_data_software_version']
        
        #self.filename_level1a = global_attrs_level1a['filename']
        self.raw_file_warning = global_attrs_level1a['raw_file_warning']
        self.labview_logfile_warning = global_attrs_level1a['labview_logfile_warning']

        if self.multiday:
            print('No meta data updates for multi-days reading (only first day is saved) !')
            for i in range(1,len(self.dates)):
                for s in self.spectrometers:
                    try:
                        calibrated_data, calibration_flags, meteo_data, global_attrs_level1a = GROSOM_library.read_level1(self.filename_level1a[s][i])
                    except FileNotFoundError:
                        print('No file for this day, skipping ', self.filename_level1a[s][i])
                    else:
                        print('Read : ', self.filename_level1a[s][i])
                        self.calibrated_data[s] = xr.concat([self.calibrated_data[s],calibrated_data], dim='time')
                        self.calibration_flags[s] = xr.concat([self.calibration_flags[s],calibration_flags], dim='time')
                        self.meteo_data[s] = xr.concat([self.meteo_data[s], meteo_data], dim='time')
                        
                    
        return self.calibrated_data, self.calibration_flags, self.meteo_data

    def clean_level1a_byFlags(self):
        '''
        cleaning the flagged level1a timestamps
        '''
        for s in self.spectrometers:
            sum_flags_da = self.calibration_flags[s].sum(dim='flags')
            good_time = sum_flags_da.where(sum_flags_da.calibration_flags == 6, drop=True).time
            self.calibrated_data[s] = self.calibrated_data[s].where(self.calibrated_data[s].time == good_time, drop=True)
            
            #self.meteo_data[s] = self.meteo_data[s].where(self.meteo_data[s].time == meteo_good_times, drop=True)

        return self.calibrated_data#, self.meteo_data

    def clean_level1a_byFlags_all(self):
        '''
        cleaning the flagged level1a timestamps
        '''
        sum_flags = dict()
        for s in self.spectrometers:
            sum_flags[s] = self.calibration_flags[s].sum(dim='flags')
        
        good_time = sum_flags['AC240'].where((sum_flags['AC240'].calibration_flags == 6)&(sum_flags['U5303'].calibration_flags == 6)&(sum_flags['USRP-A'].calibration_flags == 6),drop=True).time
        
        for s in self.spectrometers:
            self.calibrated_data[s] = self.calibrated_data[s].where(self.calibrated_data[s].time == good_time, drop=True)
            
            #self.meteo_data[s] = self.meteo_data[s].where(self.meteo_data[s].time == meteo_good_times, drop=True)

        return self.calibrated_data#, self.meteo_data
    
    def check_channel_quality(self, spectrometers):
        '''
        Based on std Tb
        see equivalent matlab function
        '''
        raise NotImplementedError()
        for s in spectrometers:
            always_bad = self.return_bad_channels(self.dates, s)

            bad_channel = self.calibrated_data[s].stdTb
    
    def find_bad_channels_stdTb(self, spectrometers, stdTb_threshold, apply_on='cal', dimension=['time','channel_idx']):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if apply_on=='cal':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.dates, spectro)
                self.calibrated_data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.calibrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.calibrated_data
        elif apply_on=='int':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.dates, spectro)
                self.integrated_data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.integrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.integrated_data
        else:
            raise ValueError()

    def add_mean_Tb(self, spectrometers, around_center = True, around_center_value = 50e6):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        # fig, ax = plt.subplots(1,1,sharex=True)
        for s in spectrometers:
            if around_center:
                if 'good_channels' in self.calibrated_data[s].data_vars:
                    clean_Tb = self.calibrated_data[s].Tb.where(self.calibrated_data[s].good_channels==1)
                    clean_Tb = clean_Tb.where(abs(self.calibrated_data[s].frequencies-self.observation_frequency)<around_center_value)
                    da_mean_Tb = clean_Tb.mean(dim='channel_idx')
                    da_median_Tb = clean_Tb.median(dim='channel_idx')
                else:
                    print('no channel filter found, be careful ! ')
                    center_Tb = self.calibrated_data[s].where(abs(self.calibrated_data[s].frequencies-self.observation_frequency)<around_center_value)
                    da_mean_Tb = center_Tb.mean(dim='channel_idx')
                    da_median_Tb = center_Tb.median(dim='channel_idx')
            else:
                if 'good_channels' in self.calibrated_data[s].data_vars:
                    clean_Tb = self.calibrated_data[s].Tb.where(self.calibrated_data[s].good_channels==1)
                    da_mean_Tb = clean_Tb.mean(dim='channel_idx')
                    da_median_Tb = clean_Tb.median(dim='channel_idx')
                else:
                    print('no channel filter found, be careful ! ')
                    da_mean_Tb = self.calibrated_data[s].Tb.mean(dim='channel_idx')
                    da_median_Tb = self.calibrated_data[s].Tb.median(dim='channel_idx')

                if any(np.abs(da_mean_Tb.data - da_median_Tb.data) > 5):
                    print('Be careful, large difference between mean and median of Tb --> big spike somewhere')

            self.calibrated_data[s] = self.calibrated_data[s].assign(mean_Tb = da_mean_Tb)
            self.calibrated_data[s] = self.calibrated_data[s].assign(median_Tb = da_median_Tb)

        #     ax.plot(self.calibrated_data[s].time, self.calibrated_data[s].mean_Tb,'.',markersize=4,label=s)
        #     ax.legend()        
        # plt.title('Mean Tb')

        return self.calibrated_data

    #def plot_variables(self, spectrometers, var)

    def add_noise_level(self, spectrometers, max_diff_level = 10):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        fig, ax = plt.subplots(1,1,sharex=True)
        for s in spectrometers:
            if 'good_channels' in self.integrated_data[s].data_vars:
                clean_Tb = self.integrated_data[s].Tb.where(self.integrated_data[s].good_channels==1)
                Tb_diff = clean_Tb.diff(dim='channel_idx')
                da_std_diff_Tb = np.std(Tb_diff.where(abs(Tb_diff)<max_diff_level),1) / np.sqrt(2)
                #da_var_diff_Tb = np.var(Tb_diff.where(abs(Tb_diff)<max_diff_level),1) / 2
            else:
                raise NotImplementedError('Identification of bad channels needed first')
                

            self.integrated_data[s] = self.integrated_data[s].assign(noise_level = da_std_diff_Tb)
            ax.plot(self.integrated_data[s].time, self.integrated_data[s].noise_level,'.', label = s)
        
            ax.legend()
        
        plt.title('Noise level')

        return self.integrated_data

    def add_binned_spectra(self, spectrometers=['AC240', 'USRP-A', 'U5303'], bin_size=[16, 80, 10], dim='chunks', df=None, plot_diff=False, corrected=False):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if df is not None:
            new_freq = np.arange(self.observation_frequency-1e9 ,self.observation_frequency+1e9, df)
            for l, s in enumerate(spectrometers):
                #new_freq = np.arange(min(self.integrated_data['U5303'].frequencies[0].data),max(self.integrated_data['U5303'].frequencies[0].data),df)
                da_binned_Tb = xr.DataArray(
                    np.ones((len(self.integrated_data[s].coords[dim]), len(new_freq)-1)),
                    dims=[dim, 'bin_freq'],
                )

                if 'good_channels' in self.integrated_data[s].data_vars:
                    clean_Tb = self.integrated_data[s].Tb.where(self.integrated_data[s].good_channels==1)
                    clean_f = self.integrated_data[s].frequencies#.where(self.integrated_data[s].good_channels==1)
                    
                    for i in range(len(self.integrated_data[s].coords[dim])):
                        digitized = np.digitize(clean_f[i], new_freq)
                        binned_f = [clean_f[i][digitized==n].mean() for n in range(1,len(new_freq))]
                        binned_Tb = [clean_Tb[i][digitized==n].mean() for n in range(1,len(new_freq))]
                        da_binned_Tb.loc[i,:] = binned_Tb
                    
                    da_binned_Tb = da_binned_Tb.assign_coords({'bin_freq':binned_f})
                else:
                    raise NotImplementedError('!')

                self.integrated_data[s] = self.integrated_data[s].assign(binned_Tb = da_binned_Tb)
        else:
            #fig, axs = plt.subplots(len(spectrometers),1,sharex=True)
            for l, s in enumerate(spectrometers):
                size = bin_size[l]
                da_binned_Tb = xr.DataArray(
                    np.ones((len(self.integrated_data[s].coords[dim]), len(self.integrated_data[s].coords['channel_idx'])//size)),
                    dims=[dim, 'bin_freq'],
                )

                if 'good_channels' in self.integrated_data[s].data_vars:
                    if corrected:
                        clean_Tb = self.integrated_data[s].Tb_corr.where(self.integrated_data[s].good_channels==1)
                    else:
                        clean_Tb = self.integrated_data[s].Tb.where(self.integrated_data[s].good_channels==1)
                    
                    clean_f = self.integrated_data[s].frequencies#.where(self.integrated_data[s].good_channels==1)

                    for i in range(len(self.integrated_data[s].coords[dim])):
                        binned_f = clean_f[i].data[-size*(clean_f[i].data.size//size):].reshape(-1,size).mean(axis=1)
                        binned_Tb = clean_Tb[i].data[-size*(clean_Tb[i].data.size//size):].reshape(-1,size).mean(axis=1) 
                        if not any(np.isnan(binned_f)):
                            saved_binned_f = binned_f

                        da_binned_Tb.loc[i,:] = binned_Tb

                    da_binned_Tb = da_binned_Tb.assign_coords({'bin_freq':saved_binned_f})

                    if plot_diff:
                        raise NotImplementedError('!')
                    #lab = str(self.integrated_data[use_basis].coords[dim].data[i])
                    #axs[l].plot(da_binned_Tb.coords['bin_freq'], Tb_diff, lw=0.2, label=lab)
                    #axs[l].set_title(s)

                    #axs[l].legend()
                else:
                    raise NotImplementedError('Identification of bad channels needed first')

                if corrected:
                    self.integrated_data[s] = self.integrated_data[s].assign(binned_Tb_corr = da_binned_Tb)
                else:
                    self.integrated_data[s] = self.integrated_data[s].assign(binned_Tb = da_binned_Tb)
            #plt.title('Tb interp diff')

        return self.integrated_data

    def add_interpolated_spectra(self, spectrometers=None, use_basis='AC240', dim='chunks', right=np.nan, left=np.nan, plot_diff=False, from_binned=True, corrected=False):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if not from_binned:
            fig, axs = plt.subplots(len(spectrometers),1,sharex=True)
            for l, s in enumerate(spectrometers):
                clean_Tb = self.integrated_data[s].Tb.where(self.integrated_data[s].good_channels==1)
                clean_f = self.integrated_data[s].frequencies.where(self.integrated_data[s].good_channels==1)
                clean_Tb_basis = self.integrated_data[use_basis].Tb.where(self.integrated_data[use_basis].good_channels==1)
                #if grid is None:
                grid = self.integrated_data[use_basis].frequencies
                #else:
                #    raise NotImplementedError('grid must be matrix --> check size')
                da_interp_Tb = xr.DataArray(
                    np.ones_like(self.integrated_data[use_basis].Tb),
                    dims=[dim, 'channels_idx'],
                )
                for i in range(len(self.integrated_data[s].Tb)):
                    interp_tb = data.interpolate(
                        grid[i].data, 
                        clean_f[i].data, 
                        clean_Tb[i].data, 
                        left, 
                        right
                        )
                    da_interp_Tb[i,:] = interp_tb
                    Tb_diff = clean_Tb_basis[i]-da_interp_Tb[i].data
                    if plot_diff:
                        lab = str(self.integrated_data[use_basis].coords[dim].data[i])
                        axs[l].plot(grid[i], Tb_diff, lw=0.2, label=lab)
                        axs[l].set_title(s)
                #da_std_diff_Tb = np.std(Tb_diff.where(abs(Tb_diff)<max_diff_level),1) / np.sqrt(2)
                #da_var_diff_Tb = np.var(Tb_diff.where(abs(Tb_diff)<max_diff_level),1) / 2
                axs[l].legend()
                plt.title('Tb interp diff')
                self.integrated_data[s] = self.integrated_data[s].assign(interpolated_Tb = da_interp_Tb)
        else:
            for l, s in enumerate(spectrometers):
                if corrected:
                    clean_Tb = self.integrated_data[s].binned_Tb_corr
                else:
                    clean_Tb = self.integrated_data[s].binned_Tb
                clean_f = self.integrated_data[s].bin_freq
                if corrected:
                    clean_Tb_basis = self.integrated_data[use_basis].binned_Tb_corr
                else:
                    clean_Tb_basis = self.integrated_data[use_basis].binned_Tb
                #if grid is None:
                grid = self.integrated_data[use_basis].bin_freq
                #else:
                #    raise NotImplementedError('grid must be matrix --> check size')
                da_interp_Tb = xr.DataArray(
                    np.ones_like(self.integrated_data[use_basis].binned_Tb),
                    dims=[dim, 'bin_freq_interp'],
                )
                for i in range(len(self.integrated_data[s].Tb)):
                    interp_tb = data.interpolate(
                        grid.data, 
                        clean_f.data, 
                        clean_Tb[i].data, 
                        left, 
                        right
                        )
                    da_interp_Tb[i,:] = interp_tb
                
                da_interp_Tb = da_interp_Tb.assign_coords({'bin_freq_interp':grid.data})
                    #Tb_diff = clean_Tb_basis[i]-da_interp_Tb[i].data
                if corrected:
                    self.integrated_data[s] = self.integrated_data[s].assign(interpolated_Tb_corr = da_interp_Tb)
                else:
                    self.integrated_data[s] = self.integrated_data[s].assign(interpolated_Tb = da_interp_Tb)

        return self.integrated_data

    def add_bias_characterization_new(self, spectrometers, dim, param_slope, around_center_value):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        ''' 

        # For each spectrometer individually
        for s in spectrometers:
            #meanTb_lc = self.integrated_data[s].Tb.where(abs(self.integrated_data[s].bin_freq_interp-self.observation_frequency)<around_center_value,drop=True).mean(dim='bin_freq_interp')
            # bias_lc = meanTb_lc-meanTb_lc_basis
            # bias_lc_fract = bias_lc/meanTb_lc_basis

            # #self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc)
            # self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc_fract)

            cleanTb = self.integrated_data[s].Tb.where(self.integrated_data[s].good_channels==1)
            # Tb_diff = cleanTb-cleanTb_basis
            # Tb_diff_fract = Tb_diff/cleanTb_basis
            right_wing_center = param_slope[s][0]
            right_wing_interval = param_slope[s][1]
            left_wing_center = param_slope[s][2]
            left_wing_interval = param_slope[s][3]

            mean_Tb_right_wing = cleanTb.where(abs(self.integrated_data[s].frequencies[0]-right_wing_center)<right_wing_interval,drop=True).mean(dim='channel_idx')
            mean_Tb_left_wing = cleanTb.where(abs(self.integrated_data[s].frequencies[0]-left_wing_center)<left_wing_interval,drop=True).mean(dim='channel_idx')

            #mean_diff_right_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-right_wing_center)<right_wing_interval,drop=True).mean(dim='bin_freq_interp')
            #mean_diff_left_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-left_wing_center)<left_wing_interval,drop=True).mean(dim='bin_freq_interp')

            slope = (mean_Tb_right_wing - mean_Tb_left_wing) / (right_wing_center - left_wing_center)

            self.integrated_data[s] = self.integrated_data[s].assign(slope_indiv = slope)
            
            self.integrated_data[s].slope_indiv.attrs['left_f']=left_wing_center
            self.integrated_data[s].slope_indiv.attrs['right_f']=right_wing_center
            
            Tb0 = mean_Tb_left_wing-slope*left_wing_center

            continuum_value_line_center = self.observation_frequency * slope + Tb0

            self.integrated_data[s] = self.integrated_data[s].assign(continuum_value_line_center = continuum_value_line_center)

            line_center_amplitude = cleanTb.where(abs(self.integrated_data[s].frequencies-self.observation_frequency)<around_center_value,drop=True).mean(dim='channel_idx')
            diff_lc_continuum = line_center_amplitude - continuum_value_line_center
            self.integrated_data[s] = self.integrated_data[s].assign(line_amplitude = diff_lc_continuum)

        return self.integrated_data

    def add_bias_characterization(self, spectrometers, use_basis, dim, param_slope, around_center_value):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        ''' 

        # Base spectro:
        cleanTb_basis = self.integrated_data[use_basis].interpolated_Tb
        meanTb_lc_basis = self.integrated_data[use_basis].interpolated_Tb.where(abs(self.integrated_data[use_basis].bin_freq_interp-self.observation_frequency)<around_center_value,drop=True).mean(dim='bin_freq_interp')
        for s in spectrometers:
            meanTb_lc = self.integrated_data[s].interpolated_Tb.where(abs(self.integrated_data[s].bin_freq_interp-self.observation_frequency)<around_center_value,drop=True).mean(dim='bin_freq_interp')
            bias_lc = meanTb_lc-meanTb_lc_basis
            bias_lc_fract = bias_lc/meanTb_lc_basis

            #self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc)
            self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc_fract)

            cleanTb = self.integrated_data[s].interpolated_Tb
            Tb_diff = cleanTb-cleanTb_basis
            Tb_diff_fract = Tb_diff/cleanTb_basis
            right_wing_center = param_slope[s][0]
            right_wing_interval = param_slope[s][1]
            left_wing_center = param_slope[s][2]
            left_wing_interval = param_slope[s][3]

            mean_diff_right_wing = Tb_diff.where(abs(Tb_diff.bin_freq_interp-right_wing_center)<right_wing_interval,drop=True).mean(dim='bin_freq_interp')
            mean_diff_left_wing = Tb_diff.where(abs(Tb_diff.bin_freq_interp-left_wing_center)<left_wing_interval,drop=True).mean(dim='bin_freq_interp')

            #mean_diff_right_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-right_wing_center)<right_wing_interval,drop=True).mean(dim='bin_freq_interp')
            #mean_diff_left_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-left_wing_center)<left_wing_interval,drop=True).mean(dim='bin_freq_interp')

            slope = (mean_diff_right_wing - mean_diff_left_wing) / (right_wing_center - left_wing_center)

            self.integrated_data[s] = self.integrated_data[s].assign(slope = slope)

            deltaTb0 = mean_diff_left_wing-slope*left_wing_center

            f_deltaTb0 = -deltaTb0 / slope

            self.integrated_data[s] = self.integrated_data[s].assign(f_deltaTb0 = f_deltaTb0)

        return self.integrated_data
    
    def add_bias_characterization_corrected(self, spectrometers, use_basis, dim, param_slope, around_center_value):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        ''' 

        # Base spectro:
        cleanTb_basis = self.integrated_data[use_basis].interpolated_Tb_corr
        meanTb_lc_basis = self.integrated_data[use_basis].interpolated_Tb_corr.where(abs(self.integrated_data[use_basis].bin_freq_interp-self.observation_frequency)<around_center_value,drop=True).mean(dim='bin_freq_interp')
        for s in spectrometers:
            meanTb_lc = self.integrated_data[s].interpolated_Tb_corr.where(abs(self.integrated_data[s].bin_freq_interp-self.observation_frequency)<around_center_value,drop=True).mean(dim='bin_freq_interp')
            bias_lc = meanTb_lc-meanTb_lc_basis
            bias_lc_fract = bias_lc/meanTb_lc_basis

            #self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc_corr = bias_lc)
            self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc_corr = bias_lc_fract)

            cleanTb = self.integrated_data[s].interpolated_Tb_corr
            Tb_diff = cleanTb-cleanTb_basis
            #Tb_diff_fract = Tb_diff/cleanTb_basis
            right_wing_center = param_slope[s][0]
            right_wing_interval = param_slope[s][1]
            left_wing_center = param_slope[s][2]
            left_wing_interval = param_slope[s][3]

            mean_diff_right_wing = Tb_diff.where(abs(Tb_diff.bin_freq_interp-right_wing_center)<right_wing_interval,drop=True).mean(dim='bin_freq_interp')
            mean_diff_left_wing = Tb_diff.where(abs(Tb_diff.bin_freq_interp-left_wing_center)<left_wing_interval,drop=True).mean(dim='bin_freq_interp')

            #mean_diff_right_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-right_wing_center)<right_wing_interval,drop=True).mean(dim='bin_freq_interp')
            #mean_diff_left_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-left_wing_center)<left_wing_interval,drop=True).mean(dim='bin_freq_interp')

            slope = (mean_diff_right_wing - mean_diff_left_wing) / (right_wing_center - left_wing_center)

            self.integrated_data[s] = self.integrated_data[s].assign(slope_corr = slope)

            deltaTb0 = mean_diff_left_wing-slope*left_wing_center

            f_deltaTb0 = -deltaTb0 / slope

            self.integrated_data[s] = self.integrated_data[s].assign(f_deltaTb0_corr = f_deltaTb0)

        return self.integrated_data

    def integrate_by_mean_Tb(self, spectrometers, Tb_chunks = [50, 100, 150, 200]):
        '''
        integration based on the mean brightness temperature in a given dataset

        output an integrated dataset with dimension (time, chunks)
        '''
        print(f"Integration with Tb chunks: {Tb_chunks}")
        self.integrated_data = dict()

        # At this point this makes no sense to integrate any meteo data...
        # this is not straightforward because the time of both datasets are different.
        self.integrated_meteo = self.meteo_data
        
        for s in spectrometers:
            chunks_num_el = []
            for i in range(len(Tb_chunks)+1):
                if i==0:
                    chunks_lim = Tb_chunks[i]
                    ds = self.calibrated_data[s].where(self.calibrated_data[s].mean_Tb < chunks_lim, drop=True)
                    chunks_num_el.append(ds.dims['time'])
                    
                    #meteo = self.meteo_data[s].where(self.calibrated_data[s].mean_Tb < chunks_lim, drop=True)
                    
                    self.integrated_data[s] = ds.mean(dim='time', skipna=True)
                    #self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                elif i==len(Tb_chunks):
                    chunks_lim = Tb_chunks[i-1]
                    ds = self.calibrated_data[s].where(self.calibrated_data[s].mean_Tb > chunks_lim, drop=True)
                    chunks_num_el.append(ds.dims['time'])
                    #meteo = self.meteo_data[s].where(self.calibrated_data[s].mean_Tb > chunks_lim, drop=True)
                    self.integrated_data[s]=xr.concat([self.integrated_data[s],ds.mean(dim='time', skipna=True)], dim='chunks')    
                    #self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='chunks')         
                else:
                    ds = self.calibrated_data[s].where(self.calibrated_data[s].mean_Tb > Tb_chunks[i-1], drop=True)
                    ds = ds.where(self.calibrated_data[s].mean_Tb < Tb_chunks[i], drop=True)
                    chunks_num_el.append(ds.dims['time'])
                    #meteo = self.meteo_data[s].where(self.calibrated_data[s].mean_Tb > Tb_chunks[i-1], drop=True)
                    #meteo = meteo.where(self.calibrated_data[s].mean_Tb < Tb_chunks[i], drop=True)
                    self.integrated_data[s]=xr.concat([self.integrated_data[s],ds.mean(dim='time', skipna=True)], dim='chunks')  
                    #self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='chunks')      
            
            self.integrated_data[s] = self.integrated_data[s].assign_coords({'chunks':np.hstack((0, Tb_chunks))})
            #self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'chunks':np.hstack((0, Tb_chunks))})
            
            # New variable with the size of the chunks
            self.integrated_data[s] = self.integrated_data[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['chunks'])
                )

            new_stdTb = np.ones(shape=self.integrated_data[s].stdTb.data.shape)
            # #noise = np.ones(shape=self.integrated_data[s].mean_Tb.data.shape)
            # #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_data[s].stdTb[i].data / np.sqrt(chunks_num_el[i])
                #noise[i] = np.std(np.diff(self.integrated_data[s].Tb[i].data[self.integrated_data[s].good_channels[i]]))
                #new_mean_stdTb[i] = np.nanmean(new_stdTb[i])
            
            # self.integrated_data[s] = self.integrated_data[s].assign(
            #     stdTb = xr.DataArray(new_stdTb, dims = ['chunks','channel_idx'])
            #     )
            
            #self.integrated_data[s] = self.integrated_data[s].assign(
            #    mean_stdTb = xr.DataArray(new_mean_stdTb, dims = ['chunks'])
            #    )

            # Drop meaningless variables:
            self.integrated_data[s] = self.integrated_data[s].drop(
                ['stdTHot','stdTSys', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )

        return self.integrated_data, self.integrated_meteo
   
    def integrate_by_mean_Tb_harmonized(self, spectrometers, Tb_chunks = [50, 100, 150, 200], use_basis='U5303'):
        '''
        integration based on the mean brightness temperature in a given dataset. In the case of the
        harmonized version, we only use 1 spectrometer to decide which time stamp we integrate.

        output an integrated dataset with dimension (time, chunks)
        '''
        print(f"Integration with Tb chunks: {Tb_chunks}")
        print(use_basis, ' is taken as basis')
        self.integrated_data = dict()
        self.integrated_meteo = dict()
        # At this point this makes no sense to integrate any meteo data...
        # this is not straightforward because the time of both datasets are different.
        #self.integrated_meteo = self.meteo_data
        
        # Use the basis spectro for defining the right timestamp.
        
        for s in spectrometers:
            chunks_num_el = []
            
            for i in range(len(Tb_chunks)+1):
                if i==0:
                    chunks_lim = Tb_chunks[i]
                    good_times = self.calibrated_data[use_basis].where(self.calibrated_data[use_basis].mean_Tb < chunks_lim, drop=True)
                    #good_times = good_times.time.values.astype('datetime64[m]')
                    #ds = self.calibrated_data[s].where(self.calibrated_data[use_basis].mean_Tb < chunks_lim, drop=True)

                    ds = self.calibrated_data[s].sel(time=good_times.time, method='nearest', drop=True)
                    chunks_num_el.append(ds.dims['time'])
                    #mean_time = ds.dims['time']
                    meteo = self.meteo_data[s].sel(time=good_times.time, method='nearest', drop=True)
                    
                    self.integrated_data[s] = ds.mean(dim='time', skipna=True)
                    self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                    #self.integrated_meteo[s]['new_prec'] = meteo.precipitation.sum(dim='time', skipna=True)

                elif i==len(Tb_chunks):
                    chunks_lim = Tb_chunks[i-1]
                    #ds = self.calibrated_data[s].where(self.calibrated_data[use_basis].mean_Tb > chunks_lim, drop=True)
                    good_times = self.calibrated_data[use_basis].where(self.calibrated_data[use_basis].mean_Tb > chunks_lim, drop=True)
                    
                    ds = self.calibrated_data[s].sel(time=good_times.time, method='nearest', drop=True)
                    chunks_num_el.append(ds.dims['time'])
                    meteo = self.meteo_data[s].sel(time=good_times.time, method='nearest', drop=True)
                    self.integrated_data[s]=xr.concat([self.integrated_data[s],ds.mean(dim='time', skipna=True)], dim='chunks')    
                    self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='chunks')     
                    #self.integrated_meteo[s].new_prec = xr.concat([self.integrated_meteo[s].new_prec,meteo.precipitation.sum(dim='time', skipna=True)], dim='chunks')    
                else:
                    good_times = self.calibrated_data[use_basis].where(self.calibrated_data[use_basis].mean_Tb > Tb_chunks[i-1], drop=True)
                    good_times = good_times.where(self.calibrated_data[use_basis].mean_Tb < Tb_chunks[i], drop=True)               
                    
                    ds = self.calibrated_data[s].sel(time=good_times.time, method='nearest', drop=True)

                    #ds = self.calibrated_data[s].where(self.calibrated_data[use_basis].mean_Tb > Tb_chunks[i-1], drop=True)
                    #ds = ds.where(self.calibrated_data[use_basis].mean_Tb < Tb_chunks[i], drop=True)
                    chunks_num_el.append(ds.dims['time'])

                    meteo = self.meteo_data[s].sel(time=good_times.time, method='nearest', drop=True)
                    
                    self.integrated_data[s]=xr.concat([self.integrated_data[s],ds.mean(dim='time', skipna=True)], dim='chunks')  
                    self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='chunks')
                    #self.integrated_meteo[s].new_prec = xr.concat([self.integrated_meteo[s].new_prec,meteo.precipitation.sum(dim='time', skipna=True)], dim='chunks')

            self.integrated_data[s] = self.integrated_data[s].assign_coords({'chunks':np.hstack((0, Tb_chunks))})
            self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'chunks':np.hstack((0, Tb_chunks))})
            
            # New variable with the size of the chunks
            self.integrated_data[s] = self.integrated_data[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['chunks'])
                )

            new_stdTb = np.ones(shape=self.integrated_data[s].stdTb.data.shape)
            # #noise = np.ones(shape=self.integrated_data[s].mean_Tb.data.shape)
            # #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_data[s].stdTb[i].data / np.sqrt(chunks_num_el[i])

            # Drop meaningless variables:
            self.integrated_data[s] = self.integrated_data[s].drop(
                ['stdTHot', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )

        return self.integrated_data, self.integrated_meteo

    def integrate_classic(self, spectrometers, freq = 1):
        '''
        integration based on the mean brightness temperature in a given dataset

        output an integrated dataset with dimension (time, chunks)
        '''        
        tod = np.arange(0, 24, freq) + freq/2
        interval = np.ones_like(tod) * freq/2

        print(f"Integration around TOD: {tod}, with interval of {interval} hours.")
        print(f"This is equivalent to an integration time of {freq} hours.")
        self.integrated_data = dict()
        self.integrated_meteo = dict()
        for s in spectrometers:
            chunks_num_el = []
            for index, (t, it) in enumerate(zip(tod, interval)):
                ds = self.calibrated_data[s].where(self.calibrated_data[s].time_of_day > t-it, drop=True)
                ds = ds.where(self.calibrated_data[s].time_of_day < t+it, drop=True)
                
                chunks_num_el.append(ds.dims['time'])

                #meteo_hour = pd.DatetimeIndex(self.meteo_data[s].time.data).hour 
                #good_meteo_times = self.meteo_data[s].time[(meteo_hour >= t-it) & (meteo_hour < t+it)]

                #meteo = self.meteo_data[s].where(self.meteo_data[s].time == good_meteo_times, drop=True)
                
                meteo = self.meteo_data[s].sel(time=ds.time, method='nearest', drop=True)

                if index == 0:
                    self.integrated_data[s] = ds.mean(dim='time', skipna=True)
                    self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                else:
                    self.integrated_data[s]=xr.concat([self.integrated_data[s],ds.mean(dim='time', skipna=True)], dim='time') 
                    self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='time') 

            self.integrated_data[s] = self.integrated_data[s].assign_coords({'time':tod})
            self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'time':tod})
            
            # New variable with the size of the chunks
            self.integrated_data[s] = self.integrated_data[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['time'])
                )

            new_stdTb = np.ones(shape=self.integrated_data[s].stdTb.data.shape)
            # noise = np.ones(shape=self.integrated_data[s].mean_Tb.data.shape)
            # #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_data[s].stdTb[i].data / np.sqrt(chunks_num_el[i])
                #noise[i] = np.std(np.diff(self.integrated_data[s].Tb[i].data[self.integrated_data[s].good_channels[i]]))
                #new_mean_stdTb[i] = np.nanmean(new_stdTb[i])
            
            # self.integrated_data[s] = self.integrated_data[s].assign(
            #     stdTb = xr.DataArray(new_stdTb, dims = ['time','channel_idx'])
            #     )
            
            #self.integrated_data[s] = self.integrated_data[s].assign(
            #    mean_stdTb = xr.DataArray(new_mean_stdTb, dims = ['chunks'])
            #    )

            # Drop meaningless variables:
            self.integrated_data[s] = self.integrated_data[s].drop(
                ['stdTHot', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )
        return self.integrated_data, self.integrated_meteo

    def integrate_by_tod(self, spectrometers, tod = [9, 16], interval= [2, 1]):
        '''
        integration based on the mean brightness temperature in a given dataset

        output an integrated dataset with dimension (time, chunks)
        '''
        print(f"Integration around TOD: {tod}, with interval of {interval} hours.")
        self.integrated_data = dict()
        self.integrated_meteo = dict()
        for s in spectrometers:
            chunks_num_el = []
            for index, (t, it) in enumerate(zip(tod, interval)):
                ds = self.calibrated_data[s].where(self.calibrated_data[s].time_of_day > t-it, drop=True)
                ds = ds.where(self.calibrated_data[s].time_of_day < t+it, drop=True)
                
                chunks_num_el.append(ds.dims['time'])
                meteo_hour = pd.DatetimeIndex(self.meteo_data[s].time.data).hour 
                good_meteo_times = self.meteo_data[s].time[(meteo_hour >= t-it) & (meteo_hour < t+it)]

                meteo = self.meteo_data[s].where(self.meteo_data[s].time == good_meteo_times, drop=True)
                
                if index == 0:
                    self.integrated_data[s] = ds.mean(dim='time', skipna=True)
                    self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                else:
                    self.integrated_data[s]=xr.concat([self.integrated_data[s],ds.mean(dim='time', skipna=True)], dim='time') 
                    self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='time') 

            self.integrated_data[s] = self.integrated_data[s].assign_coords({'time':tod})
            self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'time':tod})
            
            # New variable with the size of the chunks
            self.integrated_data[s] = self.integrated_data[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['time'])
                )

            new_stdTb = np.ones(shape=self.integrated_data[s].stdTb.data.shape)
            #noise = np.ones(shape=self.integrated_data[s].mean_Tb.data.shape)
            #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_data[s].stdTb[i].data / np.sqrt(chunks_num_el[i])

                #new_mean_stdTb[i] = np.nanmean(new_stdTb[i])
            
            # self.integrated_data[s] = self.integrated_data[s].assign(
            #     stdTb = xr.DataArray(new_stdTb, dims = ['time','channel_idx'])
            #     )

            #self.integrated_data[s] = self.integrated_data[s].assign(
            #    mean_stdTb = xr.DataArray(new_mean_stdTb, dims = ['chunks'])
            #    )

            # Drop meaningless variables:
            self.integrated_data[s] = self.integrated_data[s].drop(
                ['stdTHot','stdTSys', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )

        return self.integrated_data, self.integrated_meteo

    def integrate_by_tod_harmonized(self, spectrometers, tod = [9, 16], interval= [2, 1], use_basis='U5303'):
        '''
        integration based on the mean brightness temperature in a given dataset

        output an integrated dataset with dimension (time, chunks)
        '''
        print(f"Integration around TOD: {tod}, with interval of {interval} hours.")
        print(use_basis, ' is taken as basis for defining the integration')
        self.integrated_data = dict()
        self.integrated_meteo = dict()
        for s in spectrometers:
            chunks_num_el = []
            for index, (t, it) in enumerate(zip(tod, interval)):

                good_times = self.calibrated_data[use_basis].where(self.calibrated_data[use_basis].time_of_day > t-it, drop=True)
                good_times = good_times.where(self.calibrated_data[use_basis].time_of_day < t+it, drop=True)

                ds = self.calibrated_data[s].sel(time=good_times.time, method='nearest', drop=True)
                
                chunks_num_el.append(ds.dims['time'])

                #meteo_hour = pd.DatetimeIndex(self.meteo_data[s].time.data).hour 
                #good_meteo_times = self.meteo_data[s].time[(meteo_hour >= t-it) & (meteo_hour < t+it)]

                meteo = self.meteo_data[s].sel(time=good_times.time, method='nearest', drop=True)
                
                if index == 0:
                    self.integrated_data[s] = ds.mean(dim='time', skipna=True)
                    self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                else:
                    self.integrated_data[s]=xr.concat([self.integrated_data[s],ds.mean(dim='time', skipna=True)], dim='time') 
                    self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='time') 

            self.integrated_data[s] = self.integrated_data[s].assign_coords({'time':tod})
            self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'time':tod})
            
            # New variable with the size of the chunks
            self.integrated_data[s] = self.integrated_data[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['time'])
                )

            new_stdTb = np.ones(shape=self.integrated_data[s].stdTb.data.shape)
            #noise = np.ones(shape=self.integrated_data[s].mean_Tb.data.shape)
            #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_data[s].stdTb[i].data / np.sqrt(chunks_num_el[i])

            # Drop meaningless variables:
            self.integrated_data[s] = self.integrated_data[s].drop(
                ['stdTHot', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )

        return self.integrated_data, self.integrated_meteo

    def integrate(self, spectrometers, strategy=None, **kwargs):
        '''
        Generic function for the integration of Integration class
        '''
        if strategy is None:
            strategy = self.integration_strategy

        if strategy == 'classic':
            print('Performing simple integration')
            if self.multiday:
                raise NotImplementedError('for multiple days, this reduces to the TOD integration')
            
            return self.integrate_classic(self.spectrometers, freq=kwargs['freq'])
        elif strategy == 'TOD':
            print('Performing integration based on the time of day')
            if len(kwargs['tod']) == len(kwargs['interval']):
                
                return self.integrate_by_tod(self.spectrometers, tod=kwargs['tod'], interval=kwargs['interval'])
            else:
                raise ValueError('Number of interval must equal the number of TOD chosen')

        elif strategy == 'TOD_harmo':
            print('TOD integration harmonized')
            if len(kwargs['tod']) == len(kwargs['interval']):
                
                return self.integrate_by_tod_harmonized(self.spectrometers, tod=kwargs['tod'], interval=kwargs['interval'], use_basis=kwargs['spectro_basis'])
            else:
                raise ValueError('Number of interval must equal the number of TOD chosen')   

        elif strategy == 'meanTb':
            #print(kwargs['Tb_chunks'])
            
            return self.integrate_by_mean_Tb(self.spectrometers, Tb_chunks=kwargs['Tb_chunks'])
        
        elif strategy == 'meanTb_harmo':
            return self.integrate_by_mean_Tb_harmonized(self.spectrometers, Tb_chunks=kwargs['Tb_chunks'], use_basis=kwargs['spectro_basis'])

        elif strategy == '...':
            print('Performing integration based on ...')
            raise NotImplementedError()
        else:
            raise NotImplementedError('unkown integration strategy')

    def return_bad_channels(self, date, spectro):
        '''
        it is overwritten by each instrument classes
        '''
        return []

    def save_dataset_level1b(self, spectrometers, datasets, groups=['spectrometer1'], extra=''):
        '''
        save as netCDF level1b file
        '''
        for s in spectrometers:
            filename = self.filename_level1b[s] + extra + '.nc'
            
            for i, group_name in enumerate(groups):
                ds = datasets[i][s]
                if i==0:
                    ds.to_netcdf(
                        path=filename,
                        mode='w',
                        group=group_name,
                        format = 'NETCDF4'
                        )
                else:
                    ds.to_netcdf(
                        path=filename,
                        mode='a',
                        group=group_name,
                        format = 'NETCDF4'
                        )

            print('Successfully saved in ', filename)

    def correct_troposphere(self):
        raise NotImplementedError('Use specific correction function for the instrument')
    
    def compare_Tb_chunks(self, dim='time', idx=[0], save = False, Tb_corr = False):
        figures = list()

        for i in idx:
            figures.append(GROSOM_library.plot_Tb_chunks(self, self.integrated_data, i)) 

        if Tb_corr:
            for i in idx:
                figures.append(GROSOM_library.plot_Tb_corr_chunks(self, self.integrated_data, i)) 
        
        if save:
            raise NotImplementedError()
    
class DataRetrieval(ABC):
    def __init__(
        self,
        instrument_name=None,
        observation_frequency=None,
        spectrometers=None,
        integration_strategy='classic',
        integration_time=None,
        date=None,
        level1_folder=None,
        level2_folder=None
        ):
        
        self.instrument_name = instrument_name
        self.date = date
        #self.datestr = date.strftime('%Y_%m_%d')
        self.observation_frequency = observation_frequency
        self.spectrometers = spectrometers
        self.integration_strategy = integration_strategy
        self.int_time = integration_time
        self.level1_folder = level1_folder
        self.level2_folder = level2_folder

        # must be false for Integration
        self.multiday = False
        try:
            if len(self.date) > 1:
                self.multiday = True
        except TypeError:
            self.multiday = False
            
    
        self.calibrated_data = dict()
        self.meteo_data = dict()
        self.calibration_flags = dict()
        self.filename_level1a = dict()
        self.filename_level1b = dict()
        extra_base =''
        for s in self.spectrometers:
            self.filename_level1a[s] = list()
            self.filename_level1b[s] = list()
            if self.multiday:
                for d in self.date:
                    self.datestr = d.strftime('%Y_%m_%d')
                    self.filename_level1a[s].append(
                        os.path.join(
                    self.level1_folder,
                    self.instrument_name + "_level1a_" +
                    s + "_" + self.datestr
                    ))
                    if self.integration_strategy == 'classic':
                        if self.int_time == 1:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_" +
                            s + "_" + self.datestr
                            )
                        else:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                            s + "_" + self.datestr
                            )
                    else:
                        self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                            s + "_" + self.datestr
                            )
            else:
                self.datestr = self.date[0].strftime('%Y_%m_%d')
                if self.integration_strategy == 'classic':
                    if self.int_time == 1:
                        self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_" +
                        s + "_" + self.datestr
                        ))
                    else:
                        self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                        s + "_" + self.datestr
                        ))
                else:
                    self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                        s + "_" + self.datestr+extra_base
                        ))

    def get_hot_load_temperature(self, time):
        ''' Get hot load temperature for a specific time.
        '''
        raise NotImplementedError("abstract base class")

    def read_level1b(self, no_flag=False, meta_data=True, extra_base=None):
        ''' 
        Reading level1b dataset and completing the information on the instrument
        
        '''
        self.integrated_data = dict()
        self.flags = dict()
        self.integrated_meteo = dict()
        #self.filename_level1a = dict()
        #self.filename_level1b = dict()
        self.filename_level2 = dict()
        for i, s in enumerate(self.spectrometers):
            self.filename_level1a[s] = os.path.join(
            self.level1_folder,
            self.instrument_name + "_level1a_" +
            s + "_" + self.datestr
            )
            if self.int_time == 1:
                # self.filename_level1b[s] = os.path.join(
                # self.level1_folder,
                # self.instrument_name + "_level1b_" +
                # s + "_" + self.datestr
                # )

                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                s + "_" + self.datestr
                )
            else:
                # self.filename_level1b[s] = os.path.join(
                # self.level1_folder,
                # self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                # s + "_" + self.datestr
                # )

                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" + str(self.int_time) +"h_" +
                s + "_" + self.datestr
                )
        
            #print('reading : ', self.filename_level1b[s])
            #self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = GROSOM_library.read_level1(self._filename_lvl1b)
            try:
                self.integrated_data[s], self.flags[s], self.integrated_meteo[s], global_attrs_level1b = GROSOM_library.read_level1(self.filename_level1b[s][0], no_flag=no_flag)
                print('Read : ', self.filename_level1b[s][0])
            except FileNotFoundError:
                print('No file for this day, skipping ', self.filename_level1b[s][0])


        if meta_data:       
            # Meta data
            self.institution = global_attrs_level1b['institution']
            self.instrument_name_from_level1b = global_attrs_level1b['instrument']
            self.location = global_attrs_level1b['location']

            self.number_of_spectrometer = global_attrs_level1b['number_of_spectrometer']
            self.calibration_version = global_attrs_level1b['calibration_version']

            #self.raw_data_filename = global_attrs_level1b['raw_data_filename']
            self.raw_data_software_version = global_attrs_level1b['raw_data_software_version']

            self.filename_level1a = global_attrs_level1b['filename_level1a']
            self.raw_file_warning = global_attrs_level1b['raw_file_warning']
            self.labview_logfile_warning = global_attrs_level1b['labview_logfile_warning']

        # some information from the ds
        #self.number_of_channels = len(self.level1b_ds.channel_idx.values)

        #self.frequencies = self.level1b_ds.frequencies.values
        #self.IF = self.level1b_ds.intermediate_freq.values
        if self.multiday:
            print('No meta data updates for multi-days reading (only first day is saved) !')
            for i in range(1,len(self.date)):
                for s in self.spectrometers:
                    try:
                        calibrated_data, calibration_flags, meteo_data, global_attrs_level1b = GROSOM_library.read_level1(self.filename_level1b[s][i], no_flag=no_flag)
                    except FileNotFoundError:
                        print('No file for this day, skipping ', self.filename_level1b[s][i])
                    else:
                        print('Read : ', self.filename_level1b[s][i])
                        self.calibrated_data[s] = xr.concat([self.calibrated_data[s],calibrated_data], dim='time')
                        self.calibration_flags[s] = xr.concat([self.calibration_flags[s],calibration_flags], dim='time')
                        self.meteo_data[s] = xr.concat([self.meteo_data[s], meteo_data], dim='time')
                        
        #self.time = self.level1b_ds.time.values
        #self.number_of_time_records = len(self.time)
        return self.integrated_data, self.flags, self.integrated_meteo
        #return self.level1b_ds, self.flags, self.meteo_ds

    def read_level2(self, spectrometers=None, extra_base=''):
        ''' 
        Reading level2 dataset and completing the information on the instrument
        
        '''
        if spectrometers is None:
            spectrometers = self.spectrometers

        if self.multiday:
            self.integrated_data = dict()
            self.filename_level2 = dict()
            self.level2_data=dict()
            print('No meta data updates for multi-days reading (only first day is saved) !')
            counter = 0
            for d in self.date:
                for s in spectrometers:
                    if self.int_time == 1:
                        self.filename_level2[s] = os.path.join(
                        self.level2_folder,
                        self.instrument_name + "_level2_" +
                        s + "_" + d.strftime('%Y_%m_%d') + extra_base
                        )
                    else:
                        self.filename_level2[s] = os.path.join(
                        self.level2_folder,
                        self.instrument_name + "_level2_" + str(self.int_time)+'h_' +
                        s + "_" + d.strftime('%Y_%m_%d') + extra_base
                        )
                    print(self.filename_level2[s]) 
                    try:
                        level2_data = xr.open_dataset(
                            self.filename_level2[s] + ".nc",
                        )                    
                    except FileNotFoundError:
                        print('No file for this day, skipping ', self.filename_level2[s])
                       # level2_data[s] = xr.Dataset()
                    else:
                        if counter == 0:
                            self.level2_data[s] = level2_data
                        else:
                            self.level2_data[s] = xr.concat([self.level2_data[s],level2_data], dim='time')
                        print('Read : ', self.filename_level2[s])
                        counter = counter + 1
        else:
            self.integrated_data = dict()
            self.filename_level2 = dict()
            self.level2_data=dict()
            for s in spectrometers:
                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                s + "_" + self.datestr + extra_base
                )

                try:
                    self.level2_data[s] = xr.open_dataset(
                        self.filename_level2[s] + ".nc",
                        )
                    print('Read : ', self.filename_level2[s])
                except FileNotFoundError:
                    print('No file for this day, skipping ', self.filename_level2[s])
        return self.level2_data

    def add_bias_characterization_with_non_linearities(self, spectrometers, dim, param_slope, around_center_value):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        ''' 

        # For each spectrometer individually
        for s in spectrometers:
            theoretical_nonlinearities =  np.polyfit([80,186,292],[0,-0.20,0],deg=2)
            fitted_poly_theoretical_nonlinearities = np.poly1d(theoretical_nonlinearities) 
        
            non_lin=fitted_poly_theoretical_nonlinearities(self.integrated_data[s].Tb.data)
            
            #meanTb_lc = self.integrated_data[s].Tb.where(abs(self.integrated_data[s].bin_freq_interp-self.observation_frequency)<around_center_value,drop=True).mean(dim='bin_freq_interp')
            # bias_lc = meanTb_lc-meanTb_lc_basis
            # bias_lc_fract = bias_lc/meanTb_lc_basis

            # #self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc)
            # self.integrated_data[s] = self.integrated_data[s].assign(bias_Tb_lc = bias_lc_fract)

            cleanTb = (self.integrated_data[s].Tb+fitted_poly_theoretical_nonlinearities(self.integrated_data[s].Tb)).where(self.integrated_data[s].good_channels==1)
            # Tb_diff = cleanTb-cleanTb_basis
            # Tb_diff_fract = Tb_diff/cleanTb_basis
            right_wing_center = param_slope[s][0]
            right_wing_interval = param_slope[s][1]
            left_wing_center = param_slope[s][2]
            left_wing_interval = param_slope[s][3]

            mean_Tb_right_wing = cleanTb.where(abs(self.integrated_data[s].frequencies[0]-right_wing_center)<right_wing_interval,drop=True).mean(dim='channel_idx')
            mean_Tb_left_wing = cleanTb.where(abs(self.integrated_data[s].frequencies[0]-left_wing_center)<left_wing_interval,drop=True).mean(dim='channel_idx')

            #mean_diff_right_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-right_wing_center)<right_wing_interval,drop=True).mean(dim='bin_freq_interp')
            #mean_diff_left_wing = Tb_diff_fract.where(abs(Tb_diff_fract.bin_freq_interp-left_wing_center)<left_wing_interval,drop=True).mean(dim='bin_freq_interp')

            slope = (mean_Tb_right_wing - mean_Tb_left_wing) / (right_wing_center - left_wing_center)

            self.integrated_data[s] = self.integrated_data[s].assign(slope_indiv_nonlin = slope)
            
            self.integrated_data[s].slope_indiv.attrs['left_f']=left_wing_center
            self.integrated_data[s].slope_indiv.attrs['right_f']=right_wing_center
            
            Tb0 = mean_Tb_left_wing-slope*left_wing_center

            continuum_value_line_center = self.observation_frequency * slope + Tb0

            self.integrated_data[s] = self.integrated_data[s].assign(continuum_value_line_center_nonlin = continuum_value_line_center)

            line_center_amplitude = cleanTb.where(abs(self.integrated_data[s].frequencies-self.observation_frequency)<around_center_value,drop=True).mean(dim='channel_idx')
            diff_lc_continuum = line_center_amplitude - continuum_value_line_center
            self.integrated_data[s] = self.integrated_data[s].assign(line_amplitude_nonlin = diff_lc_continuum)

        return self.integrated_data

    def plot_level1b_TB_all(self, title='', save=False, save_name='int_spectra', idx=None):
        figures = list()
        
        if idx is None:
            figures.append(GROSOM_library.plot_Tb_all(self, self.integrated_data, title=title)) 
        else:
            figures.append(GROSOM_library.plot_Tb_selected(self, self.integrated_data, title=title, idx=idx)) 

        if save:
            save_single_pdf(self.level1_folder+'/'+save_name+self.datestr+'.pdf', figures)
            #save_pngs(self.level1_folder+'time_series_'+self.datestr+'_', figures)
    

        if save:
            save_single_pdf(self.level1_folder+'/'+save_name+self.datestr+'.pdf', figures)
            #save_pngs(self.level1_folder+'time_series_'+self.datestr+'_', figures)

    def plot_level1b_TB(self, level1b_dataset, calibration_cycle):
        plt.plot(level1b_dataset.frequencies,level1b_dataset.Tb_trop_corr[calibration_cycle])
        plt.ylim((0,200))
        pass
    
    def plot_meteo_ds_level1b_dataset(self):
        GROSOM_library.plot_meteo_level1b(self.meteo_ds)
        pass 
    
    def forward_model(self, retrieval_param):
        ''' 
        Performing a FM according to the retrieval_param dictionnary as input
        parameters
        
        '''
        return retrieval_module.forward_model(retrieval_param)
    
    def retrieve_cycle(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None, ac=None):
        ''' 
        Performing single retrieval for a given calibration cycle (defined in retrieval_param) 
        '''
        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False

        retrieval_param['ref_elevation_angle'] =  90

        return retrieval_module.retrieve_cycle(self, spectro_dataset, retrieval_param, ac_FM=ac)
    
    def retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None, ac=None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False
 
        retrieval_param['ref_elevation_angle'] = 90

        return retrieval_module.retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, ac_FM=ac)
    
    def retrieve_cycle_tropospheric_corrected_pyarts(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False
 
        retrieval_param['ref_elevation_angle'] = 90
        import retrieval_module_pyarts
        return retrieval_module_pyarts.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param)
    def test_retrieval(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        retrieval_param["binned_ch"] = False
 
        retrieval_param['ref_elevation_angle'] = 90

        return retrieval_module.test_retrieval(spectro_dataset, retrieval_param)

    def smooth_and_apply_correction(self, level1b_dataset, meteo_ds):   
        '''
        doing quick tropospheric correction as it was not saved in lvl1b

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        meteo_ds : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        return GROSOM_library.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    
    def create_binning(self, freq, tb, retrieval_param):   
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        meteo_ds : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        n_f = retrieval_param["number_of_freq_points"]
        bw_extra = 1.2*self.bandwidth

        return GROSOM_library.create_bin_vector(self.observation_frequency, freq, tb, n_f, bw_extra)

    def bin_spectrum(self, freq, tb, bin_vect):   
        '''
        

        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        meteo_ds : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''

        return GROSOM_library.bin_spectrum(freq, tb, bin_vect)

    def smooth_corr_spectra(self, level1b_dataset, retrieval_param):
        '''
        
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.
        retrieval_param : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        return GROSOM_library.smooth_corr_spectra(level1b_dataset, retrieval_param)
    
    def find_bad_channels(self, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        return GROSOM_library.find_bad_channels(self.level1b_ds, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
    
    def find_bad_channels_stdTb(self, spectrometers, stdTb_threshold, apply_on='cal', dimension=['time','channel_idx']):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if apply_on=='cal':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.date, spectro)
                self.calibrated_data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.calibrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.calibrated_data
        elif apply_on=='int':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.date, spectro)
                self.integrated_data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.integrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.integrated_data
        else:
            raise ValueError()

    def find_bad_channels_stdTb_old(self, spectro, stdTb_threshold):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        bad_channels = self.return_bad_channels(date, spectro)
        self.data[spectro] = GROSOM_library.find_bad_channels_stdTb(self.data[spectro], bad_channels, stdTb_threshold)

        return self.data[spectro]
    
    def plot_level2_from_tropospheric_corrected_spectra(self, ac, spectro_dataset, retrieval_param, title, figure_list):
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
        return GROSOM_library.plot_level2_from_tropospheric_corrected_spectra(spectro_dataset, ac, retrieval_param, title, figure_list)
    
    def plot_level2(self, ac, spectro_dataset, retrieval_param, title, figure_list):
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
        return GROSOM_library.plot_level2(spectro_dataset, ac, retrieval_param, title, figure_list)

    def plot_level2_test_retrieval(self, ac, retrieval_param, title):
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
        return GROSOM_library.plot_level2_test_retrieval(ac, retrieval_param, title)
    
    def return_bad_channels(self, date, spectro):
        '''
        it is overwritten by each instrument classes
        '''
        return []

def run_retrieval(instrument,retrieval_param):
    a=2
    '''
    In this function we call the retrieval process step-by-step
    
    instrument is actually a concrete implementation of the the abstract 
    DataRetrieval class --> will be initiated here in the future
    
    In the end, we will store this function elsewhere along with the classes definition.
    '''
    return a

# if __name__=="__main__":
    
#     instrument_name = "GROMOS"
#     date = datetime.date(2019,2,1)
#     int_time = 1

#     basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
#     basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    
#     #basename_lvl1 = "/scratch/GROSOM/Level1/"
#     #basename_lvl2 = "/scratch/GROSOM/Level2/"
    
#     line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
#     #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"

#     if instrument_name=="GROMOS":
#         from gromos_classes import IntegrationGROMOS, GROMOS_LvL2
#         calibration = IntegrationGROMOS(date, basename_lvl1, basename_lvl2, int_time)
#         instrument = GROMOS_LvL2(date, basename_lvl1, basename_lvl2, int_time)
#     elif instrument_name=="SOMORA":
#         from somora_classes import SOMORA_LvL2
#         instrument = SOMORA_LvL2(date, basename_lvl1, basename_lvl2, int_time)
#     elif instrument_name=="mopi5":
#         from mopi5_classes import MOPI5_LvL2
#         #basename_lvl1 = "/scratch/MOPI5/Level1/"
#         #basename_lvl2 = "/scratch/MOPI5/Level2/"
#         basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
#         basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
#         instrument = MOPI5_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    
#     # Dictionnary containing all EXTERNAL retrieval parameters 
#     retrieval_param = dict()

#     # type of retrieval to do:
#     # 1. tropospheric corrected
#     # 2. with h20
#     # 3. test retrieving the FM
#     retrieval_param["retrieval_type"] = 1

#     retrieval_param["obs_freq"] = instrument.observation_frequency
    
#     retrieval_param["plot_meteo_ds"] = True
#     retrieval_param["number_of_freq_points"] = 601
#     retrieval_param["irregularity_f_grid"] = 10
#     retrieval_param["show_f_grid"] = True

#     retrieval_param["z_top_sim_grid"] = 97e3
#     retrieval_param["z_bottom_sim_grid"] = 800
#     retrieval_param["z_resolution_sim_grid"] = 1e3

#     retrieval_param["z_top_ret_grid"] = 95e3
#     retrieval_param["z_bottom_ret_grid"] = 800
#     retrieval_param["z_resolution_ret_grid"] = 3e3

#     retrieval_param["z_top_ret_grid_h2o"] = 50e3
#     retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

#     retrieval_param['increased_var_factor'] = 1
#     retrieval_param['unit_var_y']  = 3


#     retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
#     retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
#     retrieval_param["apriori_O3_cov"] = 1.5e-6

#     retrieval_param['water_vapor_model'] = "H2O-PWR98"
#     #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
#     #retrieval_param["azimuth_angle"]=32
    
#     #retrieval_param['obs_freq'] = 1.4217504e11
#     retrieval_param['line_file'] = line_file
#     retrieval_param['atm'] ='fascod'
#     retrieval_param['ecmwf_store_location'] ='/scratch/ECMWF'
#     retrieval_param['extra_time_ecmwf'] = 6

#     retrieval_param['cira86_path'] = os.path.join(
#         ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')
#     # Check the structure of the file and maybe use it ?
#     #print(netCDF4.Dataset(filename+".nc").groups.keys())
    
#     # level1b_dataset, meteo_ds, global_attr_lvl1b_ds = run_retrieval(instrument,retrievalParam)
#     calibrated_data, cal_flags, cal_meteo = calibration.read_level1a()
#     data, flags, meteo = instrument.read_level1b()
#     assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'

#     if instrument_name == 'mopi5':
#         instrument = instrument.correction_function_mopi5('U5303', 290)
#         instrument.plot_comparison_mopi5_spectrometers(calibration_cycle=[0,1,2,3])
         
#     #for i,s in enumerate(instrument.spectrometer):
#     spectro = 'AC240'
#     spectro_dataset = instrument.data[spectro]
#     #retrieval_param = {**global_attrs_level1b, **retrieval_param}
#     #else :
#     #    raise ValueError('incoherent instrument definition')
    
#     #level1b_dataset = instrument.smooth_and_apply_correction(level1b_dataset, meteo_ds)
#     #instrument.return_bad_channels(date,'U5303')
#     #spectro_dataset = instrument.find_bad_channels_stdTb(spectro, 20)

#     # instrument.data = instrument.data.find_bad_channels(
#     #     spectrometer=spectro,
#     #     Tb_min=0,
#     #     Tb_max=260,
#     #     boxcar_size=128,
#     #     boxcar_thresh=7
#     # )
    
#     # instrument.plot_meteo_ds_level1b_dataset()

#     #level1b_dataset = instrument.smooth_corr_spectra(level1b_dataset, retrieval_param)
#     #f_sim, y_sim = instrument.forward_model(retrieval_param)
#     #plt.plot(f_sim, y_sim[0], level1b_dataset.frequencies.values, level1b_dataset.Tb[1].values)
    
#     retrieval_param["integration_cycle"] = 0

#     #bin_vector = instrument.create_binning(
#     #    instrument.frequencies,
#     #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
#     #    retrieval_param
#     #)

#     #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)

    
#     if retrieval_param["retrieval_type"] == 1:
#         retrieval_param["surface_altitude"] = 10e3
#         retrieval_param["observation_altitude"] =  15e3
#         ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
#         figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(ac, spectro_dataset, retrieval_param, title = 'retrieval_trop_corr')
#         level2 = ac.get_level2_xarray()
#         #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
#         save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
#     elif retrieval_param["retrieval_type"] == 2:
#         retrieval_param["surface_altitude"] = 1200
#         retrieval_param["observation_altitude"] =  1200
#         ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
#         figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='retrieval_o3_h20')
#         save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_'+retrieval_param['water_vapor_model']+'.pdf', figure_list)
#     elif retrieval_param["retrieval_type"] == 3:
#         retrieval_param["surface_altitude"] = 1500
#         retrieval_param["observation_altitude"] =  1500
#         ac, retrieval_param = instrument.test_retrieval(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
#         figure_list = instrument.plot_level2_test_retrieval(ac, retrieval_param, title ='test_retrieval_o3_h2o')
#         save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'Perrin_with_h2o.pdf', figure_list)
#     else:
#         pass


# %%
