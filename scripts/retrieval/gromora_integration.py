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
from gromora_utils import save_single_pdf

from dotenv import load_dotenv

# For ARTS, we need to specify some paths
#load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

import GROMORA_library
#import retrieval_module
#import mopi5_retrievals

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
                self.calibrated_data[s], self.calibration_flags[s], self.meteo_data[s], global_attrs_level1a = GROMORA_library.read_level1(self.filename_level1a[s][0])
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
                        calibrated_data, calibration_flags, meteo_data, global_attrs_level1a = GROMORA_library.read_level1(self.filename_level1a[s][i])
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
                self.calibrated_data[spectro] = GROMORA_library.find_bad_channels_stdTb(self.calibrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.calibrated_data
        elif apply_on=='int':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.dates, spectro)
                self.integrated_data[spectro] = GROMORA_library.find_bad_channels_stdTb(self.integrated_data[spectro], bad_channels, stdTb_threshold, dimension)
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

    def add_bias_characterization_shifted(self, spectrometers, dim, param_slope, around_center_value):
        '''
        Try a bias shifted away from the line center
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

            continuum_value_line_center = 110.88e9 * slope + Tb0

            self.integrated_data[s] = self.integrated_data[s].assign(continuum_value_line_center = continuum_value_line_center)

            line_center_amplitude = cleanTb.where(abs(self.integrated_data[s].frequencies- 110.87e9)<around_center_value,drop=True).mean(dim='channel_idx')
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
            figures.append(GROMORA_library.plot_Tb_chunks(self, self.integrated_data, i)) 

        if Tb_corr:
            for i in idx:
                figures.append(GROMORA_library.plot_Tb_corr_chunks(self, self.integrated_data, i)) 
        
        if save:
            raise NotImplementedError()
