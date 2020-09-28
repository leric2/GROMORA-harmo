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
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

import GROSOM_library
import retrieval_module
import mopi5_retrievals

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

        if len(self.dates) > 1:
            self.multiday = True
        else:
            self.multiday = False

        self.calibrated_data = dict()
        self.meteo_data = dict()
        self.calibration_flags = dict()
        self.filename_level1a = dict()
        self.filename_level1b = dict()

        for s in self.spectrometers:
            self.filename_level1a[s] = list()

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

    def read_level1a(self):
        ''' 
        Reading level1a in netCDF format
        
        '''
            
        for s in self.spectrometers:
            print('reading : ', self.filename_level1a[s][0])
            #self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = GROSOM_library.read_level1(self._filename_lvl1b)
            self.calibrated_data[s], self.calibration_flags[s], self.meteo_data[s], global_attrs_level1a = GROSOM_library.read_level1(self.filename_level1a[s][0])
        
        # Meta data
        self.institution = global_attrs_level1a['institution']
        self.instrument_name_from_level1a = global_attrs_level1a['instrument']
        self.location = global_attrs_level1a['location']

        self.number_of_spectrometer = global_attrs_level1a['number_of_spectrometer']
        self.calibration_version = global_attrs_level1a['calibration_version']

        self.raw_data_filename = global_attrs_level1a['raw_data']
        self.raw_data_software_version = global_attrs_level1a['raw_data_software_version']
        
        #self.filename_level1a = global_attrs_level1a['filename']
        self.raw_file_warning = global_attrs_level1a['raw_file_warning']
        self.labview_logfile_warning = global_attrs_level1a['labview_logfile_warning']

        if self.multiday:
            for i in range(1,len(self.dates)):
                for s in self.spectrometers:
                    print('reading : ', self.filename_level1a[s][i])
                    calibrated_data, calibration_flags, meteo_data, global_attrs_level1a = GROSOM_library.read_level1(self.filename_level1a[s][i])
                    self.calibrated_data[s] = xr.concat([self.calibrated_data[s],calibrated_data], dim='time')
                    self.calibration_flags[s] = xr.concat([self.calibration_flags[s],calibration_flags], dim='time')
                    self.meteo_data[s] = xr.concat([self.meteo_data[s], meteo_data], dim='time')
                    print('No meta data updates for multi-days reading (only first day is saved) !')

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
                self.integrated_dataset[spectro] = GROSOM_library.find_bad_channels_stdTb(self.integrated_dataset[spectro], bad_channels, stdTb_threshold, dimension)
            return self.integrated_dataset
        else:
            raise ValueError()

    def add_mean_Tb(self, spectrometers):
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
        for s in spectrometers:
            if 'good_channels' in self.calibrated_data[s].data_vars:
                clean_Tb = self.calibrated_data[s].Tb.where(self.calibrated_data[s].good_channels)
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

        return self.calibrated_data

    def integrate_by_mean_Tb(self, spectrometers, Tb_chunks = [50, 100, 150, 200]):
        '''
        integration based on the mean brightness temperature in a given dataset

        output an integrated dataset with dimension (time, chunks)
        '''
        print(f"Integration with Tb chunks: {Tb_chunks}")
        self.integrated_dataset = dict()

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
                    
                    self.integrated_dataset[s] = ds.mean(dim='time', skipna=True)
                    #self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                elif i==len(Tb_chunks):
                    chunks_lim = Tb_chunks[i-1]
                    ds = self.calibrated_data[s].where(self.calibrated_data[s].mean_Tb > chunks_lim, drop=True)
                    chunks_num_el.append(ds.dims['time'])
                    #meteo = self.meteo_data[s].where(self.calibrated_data[s].mean_Tb > chunks_lim, drop=True)
                    self.integrated_dataset[s]=xr.concat([self.integrated_dataset[s],ds.mean(dim='time', skipna=True)], dim='chunks')    
                    #self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='chunks')         
                else:
                    ds = self.calibrated_data[s].where(self.calibrated_data[s].mean_Tb > Tb_chunks[i-1], drop=True)
                    ds = ds.where(self.calibrated_data[s].mean_Tb < Tb_chunks[i], drop=True)
                    chunks_num_el.append(ds.dims['time'])
                    #meteo = self.meteo_data[s].where(self.calibrated_data[s].mean_Tb > Tb_chunks[i-1], drop=True)
                    #meteo = meteo.where(self.calibrated_data[s].mean_Tb < Tb_chunks[i], drop=True)
                    self.integrated_dataset[s]=xr.concat([self.integrated_dataset[s],ds.mean(dim='time', skipna=True)], dim='chunks')  
                    #self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='chunks')      
            
            self.integrated_dataset[s] = self.integrated_dataset[s].assign_coords({'chunks':np.hstack((0, Tb_chunks))})
            #self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'chunks':np.hstack((0, Tb_chunks))})
            
            # New variable with the size of the chunks
            self.integrated_dataset[s] = self.integrated_dataset[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['chunks'])
                )

            new_stdTb = np.ones(shape=self.integrated_dataset[s].stdTb.data.shape)
            #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_dataset[s].stdTb[i].data / chunks_num_el[i]
                #new_mean_stdTb[i] = np.nanmean(new_stdTb[i])
            
            self.integrated_dataset[s] = self.integrated_dataset[s].assign(
                stdTb = xr.DataArray(new_stdTb, dims = ['chunks','channel_idx'])
                )

            #self.integrated_dataset[s] = self.integrated_dataset[s].assign(
            #    mean_stdTb = xr.DataArray(new_mean_stdTb, dims = ['chunks'])
            #    )

            # Drop meaningless variables:
            self.integrated_dataset[s] = self.integrated_dataset[s].drop(
                ['time_min','stdTHot','stdTSys', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )

        return self.integrated_dataset, self.integrated_meteo

    def integrate_classic(self, spectrometers, freq = 1):
        '''
        integration based on the mean brightness temperature in a given dataset

        output an integrated dataset with dimension (time, chunks)
        '''        
        tod = np.arange(0, 24, freq) + freq/2
        interval = np.ones_like(tod) * freq/2

        print(f"Integration around TOD: {tod}, with interval of {interval} hours.")
        print(f"This is equivalent to an integration time of {freq} hours.")
        self.integrated_dataset = dict()
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
                    self.integrated_dataset[s] = ds.mean(dim='time', skipna=True)
                    self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                else:
                    self.integrated_dataset[s]=xr.concat([self.integrated_dataset[s],ds.mean(dim='time', skipna=True)], dim='time') 
                    self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='time') 

            self.integrated_dataset[s] = self.integrated_dataset[s].assign_coords({'time':tod})
            self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'time':tod})
            
            # New variable with the size of the chunks
            self.integrated_dataset[s] = self.integrated_dataset[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['time'])
                )

            new_stdTb = np.ones(shape=self.integrated_dataset[s].stdTb.data.shape)
            #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_dataset[s].stdTb[i].data / chunks_num_el[i]
                #new_mean_stdTb[i] = np.nanmean(new_stdTb[i])
            
            self.integrated_dataset[s] = self.integrated_dataset[s].assign(
                stdTb = xr.DataArray(new_stdTb, dims = ['time','channel_idx'])
                )

            #self.integrated_dataset[s] = self.integrated_dataset[s].assign(
            #    mean_stdTb = xr.DataArray(new_mean_stdTb, dims = ['chunks'])
            #    )

            # Drop meaningless variables:
            self.integrated_dataset[s] = self.integrated_dataset[s].drop(
                ['time_min','stdTHot','stdTSys', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )

        return self.integrated_dataset, self.integrated_meteo

    def integrate_by_tod(self, spectrometers, tod = [9, 16], interval= [2, 1]):
        '''
        integration based on the mean brightness temperature in a given dataset

        output an integrated dataset with dimension (time, chunks)
        '''
        print(f"Integration around TOD: {tod}, with interval of {interval} hours.")
        self.integrated_dataset = dict()
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
                    self.integrated_dataset[s] = ds.mean(dim='time', skipna=True)
                    self.integrated_meteo[s] = meteo.mean(dim='time', skipna=True)
                else:
                    self.integrated_dataset[s]=xr.concat([self.integrated_dataset[s],ds.mean(dim='time', skipna=True)], dim='time') 
                    self.integrated_meteo[s]=xr.concat([self.integrated_meteo[s],meteo.mean(dim='time', skipna=True)], dim='time') 

            self.integrated_dataset[s] = self.integrated_dataset[s].assign_coords({'time':tod})
            self.integrated_meteo[s] = self.integrated_meteo[s].assign_coords({'time':tod})
            
            # New variable with the size of the chunks
            self.integrated_dataset[s] = self.integrated_dataset[s].assign(
                chunk_size = xr.DataArray(chunks_num_el, dims = ['time'])
                )

            new_stdTb = np.ones(shape=self.integrated_dataset[s].stdTb.data.shape)
            #new_mean_stdTb = np.zeros(len(chunks_num_el))

            for i in range(len(chunks_num_el)):
                new_stdTb[i] = self.integrated_dataset[s].stdTb[i].data / chunks_num_el[i]
                #new_mean_stdTb[i] = np.nanmean(new_stdTb[i])
            
            self.integrated_dataset[s] = self.integrated_dataset[s].assign(
                stdTb = xr.DataArray(new_stdTb, dims = ['time','channel_idx'])
                )

            #self.integrated_dataset[s] = self.integrated_dataset[s].assign(
            #    mean_stdTb = xr.DataArray(new_mean_stdTb, dims = ['chunks'])
            #    )

            # Drop meaningless variables:
            self.integrated_dataset[s] = self.integrated_dataset[s].drop(
                ['time_min','stdTHot','stdTSys', 'mean_std_Tb','stdTRoom','number_of_hot_spectra','number_of_cold_spectra','number_of_sky_spectra','good_channels']
                )

        return self.integrated_dataset, self.integrated_meteo

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
            
        elif strategy == 'meanTb':
            #print(kwargs['Tb_chunks'])
            return self.integrate_by_mean_Tb(self.spectrometers, Tb_chunks=kwargs['Tb_chunks'])  
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
            figures.append(GROSOM_library.plot_Tb_chunks(self, self.integrated_dataset, i)) 

        if Tb_corr:
            for i in idx:
                figures.append(GROSOM_library.plot_Tb_corr_chunks(self, self.integrated_dataset, i)) 
        
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
        self.datestr = date.strftime('%Y_%m_%d')
        self.observation_frequency = observation_frequency
        self.spectrometers = spectrometers
        self.integration_strategy = integration_strategy
        self.int_time = integration_time
        self.level1_folder = level1_folder
        self.level2_folder = level2_folder

        # must be false for Integration
        self.multiday = False

        self.calibrated_data = dict()
        self.meteo_data = dict()
        self.calibration_flags = dict()
        self.filename_level1a = dict()
        self.filename_level1b = dict()

        for s in self.spectrometers:
            self.filename_level1a[s] = list()

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

    def get_hot_load_temperature(self, time):
        ''' Get hot load temperature for a specific time.
        '''
        raise NotImplementedError("abstract base class")
    
    def read_level1b(self):
        ''' 
        Reading level1b dataset and completing the information on the instrument
        
        '''
        self.data = dict()
        self.flags = dict()
        self.meteo = dict()
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
        
            print('reading : ', self.filename_level1b[s])
            #self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = GROSOM_library.read_level1(self._filename_lvl1b)
            self.data[s], self.flags[s], self.meteo[s], global_attrs_level1b = GROSOM_library.read_level1(self.filename_level1b[s])
        
        # Meta data
        self.institution = global_attrs_level1b['institution']
        self.instrument_name_from_level1b = global_attrs_level1b['instrument']
        self.location = global_attrs_level1b['location']

        self.number_of_spectrometer = global_attrs_level1b['number_of_spectrometer']
        self.calibration_version = global_attrs_level1b['calibration_version']

        self.raw_data_filename = global_attrs_level1b['raw_data_filename']
        self.raw_data_software_version = global_attrs_level1b['raw_data_software_version']
        
        self.filename_level1a = global_attrs_level1b['filename_level1a']
        self.raw_file_warning = global_attrs_level1b['raw_file_warning']
        self.labview_logfile_warning = global_attrs_level1b['labview_logfile_warning']

        # some information from the ds
        #self.number_of_channels = len(self.level1b_ds.channel_idx.values)

        #self.frequencies = self.level1b_ds.frequencies.values
        #self.IF = self.level1b_ds.intermediate_freq.values

        #self.time = self.level1b_ds.time.values
        #self.number_of_time_records = len(self.time)
        return self.data, self.flags, self.meteo
        #return self.level1b_ds, self.flags, self.meteo_ds
    
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
    
    def retrieve_cycle(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None):
        ''' 
        Performing single retrieval for a given calibration cycle (defined in retrieval_param) 
        '''
        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False

        retrieval_param['ref_elevation_angle'] = 90

        return retrieval_module.retrieve_cycle(spectro_dataset, retrieval_param)
    
    def retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False
 
        retrieval_param['ref_elevation_angle'] = 90

        return retrieval_module.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param)

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

    def find_bad_channels_stdTb(self, spectro, stdTb_threshold):
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
    
    def plot_level2_from_tropospheric_corrected_spectra(self, ac, spectro_dataset, retrieval_param, title):
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
        return GROSOM_library.plot_level2_from_tropospheric_corrected(spectro_dataset, ac, retrieval_param, title)
    
    def plot_level2(self, ac, spectro_dataset, retrieval_param, title):
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
        return GROSOM_library.plot_level2(spectro_dataset, ac, retrieval_param, title)

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

if __name__=="__main__":
    
    instrument_name = "GROMOS"
    date = datetime.date(2019,2,1)
    int_time = 1

    basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    
    #basename_lvl1 = "/scratch/GROSOM/Level1/"
    #basename_lvl2 = "/scratch/GROSOM/Level2/"
    
    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
    #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"

    if instrument_name=="GROMOS":
        from gromos_classes import IntegrationGROMOS, GROMOS_LvL2
        calibration = IntegrationGROMOS(date, basename_lvl1, basename_lvl2, int_time)
        instrument = GROMOS_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    elif instrument_name=="SOMORA":
        from somora_classes import SOMORA_LvL2
        instrument = SOMORA_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    elif instrument_name=="mopi5":
        from mopi5_classes import MOPI5_LvL2
        #basename_lvl1 = "/scratch/MOPI5/Level1/"
        #basename_lvl2 = "/scratch/MOPI5/Level2/"
        basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
        basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        instrument = MOPI5_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    
    # Dictionnary containing all EXTERNAL retrieval parameters 
    retrieval_param = dict()

    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param["retrieval_type"] = 1

    retrieval_param["obs_freq"] = instrument.observation_frequency
    
    retrieval_param["plot_meteo_ds"] = True
    retrieval_param["number_of_freq_points"] = 601
    retrieval_param["irregularity_f_grid"] = 10
    retrieval_param["show_f_grid"] = True

    retrieval_param["z_top_sim_grid"] = 97e3
    retrieval_param["z_bottom_sim_grid"] = 800
    retrieval_param["z_resolution_sim_grid"] = 1e3

    retrieval_param["z_top_ret_grid"] = 95e3
    retrieval_param["z_bottom_ret_grid"] = 800
    retrieval_param["z_resolution_ret_grid"] = 3e3

    retrieval_param["z_top_ret_grid_h2o"] = 50e3
    retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

    retrieval_param['increased_var_factor'] = 1
    retrieval_param['unit_var_y']  = 3


    retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
    retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    retrieval_param["apriori_O3_cov"] = 1.5e-6

    retrieval_param['water_vapor_model'] = "H2O-PWR98"
    #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
    #retrieval_param["azimuth_angle"]=32
    
    #retrieval_param['obs_freq'] = 1.4217504e11
    retrieval_param['line_file'] = line_file
    retrieval_param['atm'] ='fascod'
    retrieval_param['ecmwf_store_location'] ='/scratch/ECMWF'
    retrieval_param['extra_time_ecmwf'] = 6

    retrieval_param['cira86_path'] = os.path.join(
        ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')
    # Check the structure of the file and maybe use it ?
    #print(netCDF4.Dataset(filename+".nc").groups.keys())
    
    # level1b_dataset, meteo_ds, global_attr_lvl1b_ds = run_retrieval(instrument,retrievalParam)
    calibrated_data, cal_flags, cal_meteo = calibration.read_level1a()
    data, flags, meteo = instrument.read_level1b()
    assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'

    if instrument_name == 'mopi5':
        instrument = instrument.correction_function_mopi5('U5303', 290)
        instrument.plot_comparison_mopi5_spectrometers(calibration_cycle=[0,1,2,3])
         
    #for i,s in enumerate(instrument.spectrometer):
    spectro = 'AC240'
    spectro_dataset = instrument.data[spectro]
    #retrieval_param = {**global_attrs_level1b, **retrieval_param}
    #else :
    #    raise ValueError('incoherent instrument definition')
    
    #level1b_dataset = instrument.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    #instrument.return_bad_channels(date,'U5303')
    #spectro_dataset = instrument.find_bad_channels_stdTb(spectro, 20)

    # instrument.data = instrument.data.find_bad_channels(
    #     spectrometer=spectro,
    #     Tb_min=0,
    #     Tb_max=260,
    #     boxcar_size=128,
    #     boxcar_thresh=7
    # )
    
    # instrument.plot_meteo_ds_level1b_dataset()

    #level1b_dataset = instrument.smooth_corr_spectra(level1b_dataset, retrieval_param)
    #f_sim, y_sim = instrument.forward_model(retrieval_param)
    #plt.plot(f_sim, y_sim[0], level1b_dataset.frequencies.values, level1b_dataset.Tb[1].values)
    
    retrieval_param["integration_cycle"] = 0

    #bin_vector = instrument.create_binning(
    #    instrument.frequencies,
    #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
    #    retrieval_param
    #)

    #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)

    
    if retrieval_param["retrieval_type"] == 1:
        retrieval_param["surface_altitude"] = 10e3
        retrieval_param["observation_altitude"] =  15e3
        ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
        figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(ac, spectro_dataset, retrieval_param, title = 'retrieval_trop_corr')
        level2 = ac.get_level2_xarray()
        #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
        save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
    elif retrieval_param["retrieval_type"] == 2:
        retrieval_param["surface_altitude"] = 1200
        retrieval_param["observation_altitude"] =  1200
        ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
        figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='retrieval_o3_h20')
        save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_'+retrieval_param['water_vapor_model']+'.pdf', figure_list)
    elif retrieval_param["retrieval_type"] == 3:
        retrieval_param["surface_altitude"] = 1500
        retrieval_param["observation_altitude"] =  1500
        ac, retrieval_param = instrument.test_retrieval(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
        figure_list = instrument.plot_level2_test_retrieval(ac, retrieval_param, title ='test_retrieval_o3_h2o')
        save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'Perrin_with_h2o.pdf', figure_list)
    else:
        pass


# %%