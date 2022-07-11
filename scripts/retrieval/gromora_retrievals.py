#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24.11.2021

@author: eric

Main class definition for GROMORA retrieval:

This file implements the abstract "DataRetrieval" class for GROMOS and SOMORA retrievals.

It is initiated by the main instrument classes (either gromos or somora class) and it contains 
all the core function to perform a retrievals from on of these 2 instruments. 
"""

from abc import ABC
import os, time, datetime 

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FuncFormatter, AutoMinorLocator, FormatStrFormatter)

import gromora_atmosphere
import GROMORA_library 

from gromora_time import get_LST_from_GROMORA, mjd2k_date
from gromora_utils import save_single_pdf, sideband_response_theory

from retrievals import arts
from retrievals import covmat

from retrievals.arts.atmosphere import p2z_simple, z2p_simple
from retrievals.data import p_interpolate

from pyarts.workspace import arts_agenda

from GROMORA_library import cmap

class DataRetrieval(ABC):
    '''
    The __init__ function is basically implementing a new class and is automatically called at each new retrieval.

    It initiates automatically the filenames of the level 1a and 1b.
    
    
    '''
    def __init__(
        self,
        instrument_name=None,
        observation_frequency=None,
        spectrometers=None,
        integration_strategy='classic',
        integration_time=None,
        date=None,
        level1_folder=None,
        level2_folder=None,
        extra_base=None
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
        #extra_base =''
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
                    s + "_" + self.datestr + extra_base
                    ))
                    if self.integration_strategy == 'classic':
                        if self.int_time == 1:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_" +
                            s + "_" + self.datestr + extra_base
                            )
                        else:
                            self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                            s + "_" + self.datestr + extra_base
                            )
                    else:
                        self.filename_level1b[s] = os.path.join(
                            self.level1_folder,
                            self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                            s + "_" + self.datestr + extra_base
                            )
            else:
                self.datestr = self.date.strftime('%Y_%m_%d')
                if self.integration_strategy == 'classic':
                    if self.int_time == 1:
                        self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_" +
                        s + "_" + self.datestr + extra_base
                        ))
                    else:
                        self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                        s + "_" + self.datestr + extra_base
                        ))
                else:
                    self.datestr = self.date.strftime('%Y_%m_%d')
                    self.filename_level1b[s].append(os.path.join(
                        self.level1_folder,
                        self.instrument_name + "_level1b_"+ self.integration_strategy + '_' +
                        s + "_" + self.datestr + extra_base
                        ))

    def correct_troposphere(self, spectrometers, dim, method='Ingold_v1'):
        '''
        Correction function for the troposphere. 
        
        Invidual correction for each spectrometers specified !
        '''
        return GROMORA_library.correct_troposphere(self, spectrometers, dim, method=method)

    #####################################################################################################################
    # Level 1b reading
    #####################################################################################################################
    def read_level1b(self, no_flag=False, meta_data=True, extra_base=None):
        """Reading level1b dataset and completing the meta-information on the instrument.


        Args:
            no_flag (bool, optional): option to avoid reading level 1 flags. Defaults to False.
            meta_data (bool, optional): option to read the meta data from level 1. Defaults to True.
            extra_base (str, optional): extra string in the level 1 filenames. Defaults to None.

        Returns:
            multiple xarray dataset: the content of the level 1b file.
        """
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
                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                s + "_" + self.datestr + extra_base
                )
            else:
                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" + str(self.int_time) +"h_" +
                s + "_" + self.datestr + extra_base
                )

            try:
                # The actual reading function
                self.integrated_data[s], self.flags[s], self.integrated_meteo[s], global_attrs_level1b = GROMORA_library.read_level1(self.filename_level1b[s][0], no_flag=no_flag)
                print('Read : ', self.filename_level1b[s][0])
            except FileNotFoundError:
                print('No file for this day, skipping ', self.filename_level1b[s][0])


        if meta_data:       
            # Meta data
            #self.source = global_attrs_level1b['source'] 
            self.name_PI = global_attrs_level1b['name'] 
            self.mail_PI = global_attrs_level1b['mail'] 
            self.contact = global_attrs_level1b['contact'] 

            #self.institution = global_attrs_level1b['institution']
            self.instrument_name_from_level1b = global_attrs_level1b['instrument']
            self.location = global_attrs_level1b['location']

            self.number_of_spectrometer = global_attrs_level1b['number_of_spectrometer']
            self.calibration_version = global_attrs_level1b['calibration_version']

            self.raw_data_filename = global_attrs_level1b['raw_filename']
            self.raw_data_software_version = global_attrs_level1b['raw_data_software_version']
            self.filename_level1b = global_attrs_level1b['filename']

            self.filename_level1a = global_attrs_level1b['filename_level1a']
            self.raw_file_warning = global_attrs_level1b['raw_file_warning']
            self.labview_logfile_warning = global_attrs_level1b['labview_logfile_warning']

            self.filtering_of_calibrated_spectra = global_attrs_level1b['filtering_of_calibrated_spectra']
            self.outlier_detection = global_attrs_level1b['outlier_detection']


        # some information from the ds
        #self.number_of_channels = len(self.level1b_ds.channel_idx.values)

        #self.frequencies = self.level1b_ds.frequencies.values
        #self.IF = self.level1b_ds.intermediate_freq.values
        if self.multiday:
            print('No meta data updates for multi-days reading (only first day is saved) !')
            for i in range(1,len(self.date)):
                for s in self.spectrometers:
                    try:
                        calibrated_data, calibration_flags, meteo_data, global_attrs_level1b = GROMORA_library.read_level1(self.filename_level1b[s][i], no_flag=no_flag)
                    except FileNotFoundError:
                        print('No file for this day, skipping ', self.filename_level1b[s][i])
                    else:
                        print('Read : ', self.filename_level1b[s][i])
                        self.calibrated_data[s] = xr.concat([self.calibrated_data[s],calibrated_data], dim='time')
                        self.calibration_flags[s] = xr.concat([self.calibration_flags[s],calibration_flags], dim='time')
                        self.meteo_data[s] = xr.concat([self.meteo_data[s], meteo_data], dim='time')
                        
        return self.integrated_data, self.flags, self.integrated_meteo

    #####################################################################################################################
    # Level 2 reading function
    #####################################################################################################################
    def read_level2(self, spectrometers=None, extra_base=''):
        """Reading level2 dataset and completing the information on the instrument


        Args:
            spectrometers (list of str, optional): list of spectrometer to read. Defaults to None: the default ones
            extra_base (str, optional): extra string to append to filename. Defaults to ''.

        Returns:
            xarray dataset: the content of the level 2
        """
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
                            decode_times=True,
                            decode_coords=True,
                            use_cftime=False,
                        )                    
                    except FileNotFoundError:
                        print('No file for this day, skipping ', self.filename_level2[s])
                       # level2_data[s] = xr.Dataset()
                    else:
                        if counter == 0:
                            self.level2_data[s] = level2_data
                        else:
                            self.level2_data[s] = xr.merge([self.level2_data[s],level2_data])
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
                        decode_times=True,
                        decode_coords=True,
                        use_cftime=False,
                        )
                    print('Read : ', self.filename_level2[s])
                except FileNotFoundError:
                    print('No file for this day, skipping ', self.filename_level2[s])
        return self.level2_data


    #####################################################################################################################
    # Definition of the retrieval_param dictionary
    #####################################################################################################################
    def define_retrieval_param(self, retrieval_param):
        """This function fills the dictionary containing all parameters for GROMOS retrievals.

        Args:
            retrieval_param (dict): a dict containing all required parameters, names, folder path etc for running a GROMORA retrievals.
        """
        retrieval_param["obs_freq"] = self.observation_frequency

        ########################################################
        # Sensor related parameter:
        retrieval_param['sensor'] = 'FFT_SB'
        retrieval_param['SB_bias'] = 0
        
        retrieval_param['sideband_response'] = 'theory'
        retrieval_param['use_all_channels'] = True
        retrieval_param['window_corrected_spectrum'] = True
        retrieval_param["f_shift"] = 0 

        # Frequency grid for the simulation:
        retrieval_param["number_of_freq_points"] = 1201
        retrieval_param["irregularity_f_grid"] = 45

        ########################################################
        # Pressure grids
        # Computed from altitudes values (all values in meters)

        # Altitude of the surface and radiometer
        retrieval_param["surface_altitude"] = 1000
        retrieval_param["observation_altitude"] =  1000   
        
        # Pressure grid for the simulation (FM)
        retrieval_param["z_top_sim_grid"] = 112e3
        retrieval_param["z_bottom_sim_grid"] = 600 #600
        retrieval_param["z_resolution_sim_grid"] = 2e3

        # Pressure grid for the retrievals 
        retrieval_param["retrieval_grid_type"] = 'altitude'
        retrieval_param["z_top_ret_grid"] = 95e3
        retrieval_param["z_bottom_ret_grid"] = 1e3
        retrieval_param["z_resolution_ret_grid"] = 2e3

        # Pressure value for the continuum retrievals
        retrieval_param["retrieval_h2o_grid_type"] = 'pressure'
        # retrieval_param["z_top_ret_grid_h2o"] = 20e3
        # retrieval_param["z_bottom_ret_grid_h2o"] = 1e3
        # retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

        retrieval_param["h2o_pressure"] = [500e2]

        # Offset for the pointing angle -> instrument_class (usually 0)
        retrieval_param['pointing_angle_corr'] = self.correct_pointing(retrieval_param)


        ########################################################
        # Species definition and spectroscopy
        retrieval_param['spectroscopy_type'] = 'XML'
        # 'H2O-PWR98' #H2O-ContMPM93 'H2O-PWR98, H2O' #"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
        retrieval_param['water_vapor_model'] = 'H2O-PWR98'
        retrieval_param['o2_model'] = 'O2-PWR93'  # 'O2-MPM93'
        # 'N2-SelfContMPM93'
        retrieval_param['n2_model'] = 'N2-SelfContStandardType'

        retrieval_param['selected_species'] = ['O3', retrieval_param['water_vapor_model'],
                                               retrieval_param['o2_model'], retrieval_param['n2_model']]

        #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
        # retrieval_param["azimuth_angle"]=32
        line_file = retrieval_param['ARTS_DATA_PATH'] +"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
        #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"
        #line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all_species.par'
        #line_file = '/home/es19m597/Documents/GROMORA/InputsRetrievals/Hitran_f_modified.xml'
        retrieval_param['line_file'] = line_file
        ########################################################
        # Atmosphere for the simulation
        retrieval_param['atm'] = 'ecmwf_cira86'  # fascod  ecmwf_cira86
        # max_diff, simple_stack_corr, simple, max_diff_surf
        retrieval_param['ptz_merge_method'] = 'max_diff'
        retrieval_param['ptz_merge_max_Tdiff'] = 5
        retrieval_param['ecmwf_store_location'] = '/storage/tub/atmosphere/ecmwf/locations/'+self.location #+str(retrieval_param['date'].year)
        #retrieval_param['ecmwf_store_location'] ='/home/eric/Documents/PhD/ECMWF'
        retrieval_param['cira86_path'] = os.path.join(retrieval_param['ARTS_DATA_PATH'], 'planets/Earth/CIRA86/monthly')
        # time interval [h] around which we collect ECMWF profile
        retrieval_param['extra_time_ecmwf'] = 3.5

        ########################################################
        # A priori and covariances
        retrieval_param['o3_apriori'] = 'waccm_monthly'#'waccm_monthly'
        # 'waccm_yearly_scaled'low_alt_ratio
        retrieval_param['h2o_apriori'] = 'ecmwf'  # 'ecmwf' # 'fascod_extended'
        retrieval_param['o3_apriori_covariance'] = 'sinefit_optimized' # 'low_alt_ratio_optimized' #low_alt_ratio_optimized
        retrieval_param['waccm_file'] = retrieval_param["GROMORA_FOLDER"] + '/files/waccm_o3_climatology.nc' #'/storage/tub/atmosphere/WACCM/waccm_o3_climatology.nc'
        retrieval_param["apriori_O3_cov"] = 1e-6  # 1e-6
        retrieval_param["apriori_H2O_stdDev"] = 1  # 6e-4 12e-4 0.5 16e-4

        ########################################################
        # Polyfit retrieval
        retrieval_param['poly_order'] = 2
        retrieval_param['covmat_polyfit_0'] = 0.1
        retrieval_param['covmat_polyfit_1'] = 0.5
        retrieval_param['covmat_polyfit_2'] = 0.5

        ########################################################
        # Type of noise covariance to use
        retrieval_param['noise_covariance']  = 'noise_level'
        # factor to increase the noise variance
        retrieval_param['increased_var_factor'] = 1

        ########################################################
        #  Baseline retrievals
        retrieval_param['sinefit_periods'] = self.baseline_period(retrieval_param) # np.array([400e6]) #np.array([319e6]) 160e6 110e6 119e6 113.5, 66e6, 45e6
        retrieval_param['sinefit_covmat'] = len(retrieval_param['sinefit_periods']) * [np.array([0.1, 0.1])] 
        retrieval_param["binned_ch"] = False

        # OEM parameters, see https://atmtools.github.io/arts-docs-2.4/docserver/methods/OEM.html 
        retrieval_param['oem_method'] = 'gn'# 'lm'
        retrieval_param['max_iter'] = 10
        retrieval_param['stop_dx'] = 20.1
        retrieval_param['lm_ga_setting']=[10, 2.0, 2.0, 100, 1, 99]

        return retrieval_param
    
    ############################################################################################################################
    # Covariances definition
    ############################################################################################################################
    def set_continuum_apriori_covariance(self, retrieval_param, p_grid_retrieval_h2o):
        """Function to built the continuum apriori covariance matrix


        Args:
            retrieval_param (dict): dictionary with all retrieval parameters 
            p_grid_retrieval_h2o (numpy array): the pressure grid retrieval for the continuum

        Returns:
            sparse covmat: covariance matrix for continuum retrievals
        """
        # sx_water = covmat.covmat_diagonal_sparse(
        #     retrieval_param["apriori_H2O_stdDev"] * np.ones_like(p_grid_retrieval_h2o))
        #sx_water = covmat.covmat_diagonal_sparse((z_grid_retrieval_h2o/z_grid_retrieval_h2o[0])*retrieval_param["apriori_H2O_stdDev"])

        # sx_o2 = covmat.covmat_diagonal_sparse(
        #     retrieval_param["apriori_o2_stdDev"] * np.ones_like(p_grid_retrieval))
        # sx_n2 = covmat.covmat_diagonal_sparse(
        #     retrieval_param["apriori_n2_stdDev"] * np.ones_like(p_grid_retrieval))
        sx_water = covmat.covmat_1d_sparse(
            grid1=np.log10(p_grid_retrieval_h2o),
            sigma1=retrieval_param["apriori_H2O_stdDev"] *
            np.ones_like(p_grid_retrieval_h2o),
            cl1=0.2*np.ones_like(p_grid_retrieval_h2o),
            fname="lin",
            cutoff=0
        )
        # std_h2o = 0.3*retrieval_param['test_apriori'][:,0]
        # plt.plot(std_h2o)
        # sx_water = covmat.covmat_1d_sparse(
        #     grid1=np.log10(p_grid_retrieval),
        #     sigma1=std_h2o,
        #     cl1=np.ones_like(p_grid_retrieval),
        #     fname="lin",
        #     cutoff=0
        # )
        #plt.matshow(sx_water.todense())
        #plt.colorbar()
        return sx_water

    def set_o3_apriori_covariance(self, retrieval_param, p_grid_retrieval):
        """Function to built the ozone apriori covariance

        Args:
            retrieval_param (dict): dictionary with all retrieval parameters 
            p_grid_retrieval (numpy array): the pressure retrieval grid

        Raises:
            ValueError: raised if the type of ozone a priori covariance saved into retrieval_param['o3_apriori_covariance'] is not recognized.

        Returns:
            sparse covmat: the a priori covariance matrix for ozone
        """
        if retrieval_param['o3_apriori_covariance']=='waccm':
            ds_waccm = gromora_atmosphere.read_waccm(retrieval_param, extra_day=1)
            sigma_o3 = p_interpolate(p_grid_retrieval, ds_waccm.p.data, ds_waccm.o3_std.data)
        elif retrieval_param['o3_apriori_covariance']=='constant':
            sigma_o3 = retrieval_param["apriori_O3_cov"]*np.ones_like(p_grid_retrieval)
        elif retrieval_param['o3_apriori_covariance']=='constant_ratio':
            sigma_o3 = retrieval_param["apriori_O3_cov"]*retrieval_param['test_apriori'][:,0]

        elif retrieval_param['o3_apriori_covariance']=='jump':
            sigma_o3 = retrieval_param["apriori_O3_cov"]*np.ones_like(p_grid_retrieval)
            sigma_o3[0:4] = 1e-7

        elif retrieval_param['o3_apriori_covariance']=='waccm_smooth_scaled':
            ds_waccm = gromora_atmosphere.read_waccm(retrieval_param, extra_day=10)
            #smoothed_std = np.convolve(ds_waccm.o3_std.data, np.ones(8)/8, mode ='same')
            smoothed_std = ds_waccm.o3_std.data
            smoothed_std[ds_waccm.p.data<100] = 0.8e-6
            smoothed_std[ds_waccm.p.data<1] = 0.4e-6
            smoothed_std = np.convolve(smoothed_std, np.ones(12)/12, mode ='same')

            sigma_o3 = p_interpolate(p_grid_retrieval, ds_waccm.p.data, smoothed_std)
            #sigma_o3 = 1e-6*sigma_o3/max(sigma_o3)

            #sigma_o3[p_grid_retrieval<100] = 0.8e-6
            #sigma_o3[p_grid_retrieval<1] = 0.2e-6
            #sigma_o3[p_grid_retrieval<p_grid_retrieval[np.where(sigma_o3 == max(sigma_o3))]] = 1e-6
            plt.semilogy(1e6*sigma_o3, 1e-2*p_grid_retrieval)
            plt.gca().invert_yaxis()
        elif retrieval_param['o3_apriori_covariance']=='waccm_monthly_scaled':
            ds_waccm = gromora_atmosphere.read_waccm_monthly(retrieval_param)
            smoothed_std = ds_waccm.o3_std.data
            smoothed_std[ds_waccm.p.data<100] = 1e-6
            smoothed_std[ds_waccm.p.data<1] = 0.4e-6
            smoothed_std = np.convolve(smoothed_std, np.ones(15)/15, mode ='same')

            sigma_o3 = p_interpolate(p_grid_retrieval, ds_waccm.p.data, smoothed_std)
            sigma_o3 = 1e-6*sigma_o3/max(sigma_o3)

        elif retrieval_param['o3_apriori_covariance']=='waccm_yearly_scaled':
            ds_waccm = gromora_atmosphere.read_waccm_yearly(retrieval_param['waccm_file'] , retrieval_param["time"])
            smoothed_std = ds_waccm.o3_std.data
            smoothed_std[ds_waccm.p.data<2000] = 0.8e-6
            smoothed_std[ds_waccm.p.data<1] = 0.4e-6
            smoothed_std = np.convolve(smoothed_std, np.ones(12)/12, mode ='same')

            sigma_o3 = p_interpolate(p_grid_retrieval, ds_waccm.p.data, smoothed_std)
            sigma_o3 = 1e-6*sigma_o3/max(sigma_o3)

            # plt.semilogy(1e6*sigma_o3, 1e-2*p_grid_retrieval)
            # plt.gca().invert_yaxis()
            # plt.semilogy(0.1*1e6*ds_waccm.o3.data, 1e-2*ds_waccm.p.data)
        elif retrieval_param['o3_apriori_covariance']=='low_alt_ratio':
            ds_waccm = gromora_atmosphere.read_waccm_yearly(retrieval_param['waccm_file'] , retrieval_param["time"])
            smoothed_std = ds_waccm.o3_std.data

            smoothed_std[ds_waccm.p.data>1000] = 0.15*ds_waccm.o3.where(ds_waccm.p>1000, drop=True).data
            smoothed_std[ds_waccm.p.data>4000] = 0.1*ds_waccm.o3.where(ds_waccm.p>4000, drop=True).data
            #smoothed_std[ds_waccm.p.data<4000] = 0.5e-6
            smoothed_std[ds_waccm.p.data<1000] = 1e-6
            smoothed_std[ds_waccm.p.data<5] = 1e-6#0.8e-6
            smoothed_std[ds_waccm.p.data<1] = 1e-6#0.6e-6
            smoothed_std = np.convolve(smoothed_std, np.ones(16)/16, mode ='same')

            sigma_o3 = p_interpolate(p_grid_retrieval, ds_waccm.p.data, smoothed_std)
            #sigma_o3 = 1e-6*sigma_o3/max(sigma_o3)

        elif retrieval_param['o3_apriori_covariance']=='low_alt_ratio_optimized':
            ds_waccm = gromora_atmosphere.read_waccm_yearly(retrieval_param['waccm_file'] , retrieval_param["time"])
            smoothed_std = ds_waccm.o3_std.data
            max_o3_p = ds_waccm.p.where(ds_waccm.o3 == max(ds_waccm.o3), drop=True).data

            smoothed_std[ds_waccm.p.data>max_o3_p] = 0.12*ds_waccm.o3.where(ds_waccm.p>max_o3_p, drop=True).data
            smoothed_std = 1e-6*smoothed_std/max(smoothed_std)
            #smoothed_std[ds_waccm.p.data>4000] = 0.1*ds_waccm.o3.where(ds_waccm.p>4000, drop=True).data
            #smoothed_std[ds_waccm.p.data<4000] = 0.5e-6
            smoothed_std[ds_waccm.p.data<=max_o3_p] = 1e-6
            smoothed_std = np.convolve(smoothed_std, np.ones(16)/16, mode ='same')

            sigma_o3 = p_interpolate(p_grid_retrieval, ds_waccm.p.data, smoothed_std)
            #sigma_o3 = 1e-6*sigma_o3/max(sigma_o3)      
        elif retrieval_param['o3_apriori_covariance']=='sinefit_optimized':
            ds_waccm = gromora_atmosphere.read_waccm_yearly(retrieval_param['waccm_file'] , retrieval_param["time"])
            smoothed_std = ds_waccm.o3_std.data
            max_o3_p = ds_waccm.p.where(ds_waccm.o3 == max(ds_waccm.o3), drop=True).data

            smoothed_std[ds_waccm.p.data>max_o3_p] = 0.2*ds_waccm.o3.where(ds_waccm.p>max_o3_p, drop=True).data
            smoothed_std = 1e-6*smoothed_std/max(smoothed_std)
            #smoothed_std[ds_waccm.p.data>4000] = 0.1*ds_waccm.o3.where(ds_waccm.p>4000, drop=True).data
            smoothed_std[ds_waccm.p.data<=max_o3_p] = 1e-6
            smoothed_std[ds_waccm.p.data<0.1] = 0.8e-6
            smoothed_std[ds_waccm.p.data<0.05] = 0.5e-6
            smoothed_std = np.convolve(smoothed_std, np.ones(16)/16, mode ='same')

            sigma_o3 = p_interpolate(p_grid_retrieval, ds_waccm.p.data, smoothed_std)
        else:
            raise ValueError('Please select another option for o3 apriori covariance matrix !')

        sx = covmat.covmat_1d_sparse(
            grid1=np.log10(p_grid_retrieval),
            sigma1= sigma_o3,
            cl1=1* np.ones_like(p_grid_retrieval),
            fname="exp",
            cutoff=0.1
        )

        if retrieval_param['plot_o3_apriori_covariance']:
            plt.figure()
            plt.semilogy(1e6*sigma_o3, p_grid_retrieval)
            plt.gca().invert_yaxis()
            # #plt.ylim((100000,100))
            # plt.semilogy(0.1*1e6*ds_waccm.o3.data, ds_waccm.p.data)

            plt.figure()
            plt.matshow(sx.todense())
            plt.colorbar()  

            #np.save('/home/es19m597/Documents/GROMORA/Data/apriori_cov', sx.todense())

        return sx

    ############################################################################################################################
    # Sensor definition
    ############################################################################################################################
    def define_sensor(self, retrieval_param, ds_freq, ds_df):
        """Definition of the GROMORA sensor type to consider for the retrieval

        Args:
            retrieval_param (dict): dictionary with all retrieval parameters 
            ds_freq (): frequency measurement vector
            ds_df (): frequency resolution of the spectrometer

        Raises:
            ValueError: type of sensors not recognized (retrieval_param['sensor'])

        Returns:
            arts.Sensor(): the sensor
        """
        # Defining our sensors
        if retrieval_param['sensor']=='FFT': 
            sensor = arts.SensorFFT(ds_freq+retrieval_param["f_shift"], ds_df)
        elif retrieval_param['sensor']=='FFT_SB': #np.concatenate((ds_freq,np.arange(148.974e9,149.977e9,1e6)))
            #intermediate_freq = 1e9*np.array([-4.101, -3.7,-3.1,3.1, 3.7,4.101])
            lo = self.lo
            if self.instrument_name == 'GROMOS':
                deltaZ = (20.04e-3 - 0.1e-3) - retrieval_param['SB_bias']
                #deltaZ1 = 20.95e-3
                lsb = 1e9*np.array([-4.101, -3.7,-3.1])
                usb= -np.flip(lsb)
                lsb_all = np.arange(-4.101e9, -3.1e9, 10e6)
                usb_all = -np.flip(lsb_all)
                plot_freq = np.arange(-4.101e9, +4.101e9, 10e6)
            elif self.instrument_name == 'SOMORA':
                deltaZ = 11.5e-3 - retrieval_param['SB_bias']
                #deltaZ1 = 20.95e-3
                lsb = 1e9*np.array([-7.601, -7.1 ,-6.6])
                usb= -np.flip(lsb)
                lsb_all = np.arange(-7.601e9, -6.6e9, 10e6)
                usb_all = -np.flip(lsb_all)
                plot_freq = np.arange(-7.601e9, +6.6e9, 10e6)
            else:
                print('SB ratio not implement for this instrument !')

            if retrieval_param['sideband_response']=='constant':
                lsb_response = np.array([1,1,1])
                usb_response = np.array([0.05,0.05,0.05])
                intermediate_freq = np.concatenate((lsb,usb))
                sideband_response = np.concatenate((lsb_response,usb_response))
            elif retrieval_param['sideband_response']=='constant_normalized':
                lsb_response = np.array([1,1,1])
                usb_response = np.array([0.05,0.05,0.05])
                lsb_response = lsb_response / (usb_response+lsb_response)
                usb_response = usb_response / (usb_response+lsb_response)
                intermediate_freq = np.concatenate((lsb,usb))
                sideband_response = np.concatenate((lsb_response,usb_response))
            elif retrieval_param['sideband_response']=='quad':
                lsb_coeff= np.polyfit(lsb, np.array([1,1,1]), deg=1)
                lsb_response = np.polyval(lsb_coeff, lsb_all)
                usb_coeff = np.polyfit(usb, np.array([0.01,0.01,0.01]), deg=1)
                usb_response = np.polyval(usb_coeff, usb_all)
                intermediate_freq = np.concatenate((lsb_all,usb_all))
                sideband_response = np.concatenate((lsb_response,usb_response))
            elif retrieval_param['sideband_response']=='theory':          
                intermediate_freq = np.concatenate((lsb_all,usb_all))
                sideband_response = sideband_response_theory(deltaZ, lo + intermediate_freq, polarisation_change=True)
            elif retrieval_param['sideband_response']=='theory_normalized':
                lower_sideband_response = sideband_response_theory(deltaZ, lo + lsb_all, polarisation_change=True)
                upper_sideband_response = sideband_response_theory(deltaZ, lo + usb_all, polarisation_change=True)
                lower_sideband_response = lower_sideband_response/(lower_sideband_response+upper_sideband_response)
                upper_sideband_response = upper_sideband_response/(lower_sideband_response+upper_sideband_response)
                intermediate_freq = np.concatenate((lsb_all,usb_all))
                sideband_response = np.concatenate((lower_sideband_response,upper_sideband_response))
            else: 
                print('select SB response type')

            #sideband_response =  np.array([1,1,1,0,0,0])
            #
            # plt.plot(1e-9*(lo+plot_freq), 100*sideband_response_theory(deltaZ, lo + plot_freq, polarisation_change=True))
            # #plt.plot(plot_freq, sideband_response_theory(deltaZ1, lo + plot_freq))
            # plt.plot(1e-9*(lo+intermediate_freq), 100*sideband_response,'x')
            # plt.xlabel('RF [GHz]')
            # plt.ylabel('MPI transmission [%]')
            # plt.grid()
            # plt.title('MPI transmission, '+retrieval_param['sideband_response']+', '+instrument.instrument_name)

            sensor = arts.SensorFFT_Sideband(ds_freq,
                ds_df, 
                num_channels=10,
                lo_freq = lo, 
                sideband_mode='lower', 
                intermediate_freq= intermediate_freq,
                sideband_response=sideband_response
            )
        elif retrieval_param['sensor']=='OFF':
            sensor = arts.SensorOff()
        else:
            raise ValueError('Please specify a valid sensor type !')
        
        return sensor

    ############################################################################################################################
    # Grids definition
    ############################################################################################################################
    # Frequency grid:
    def make_f_grid(self, retrieval_param): 
        """Function to create the frequency grid for the retrievals     

        Args:
            retrieval_param (dict): dictionary with all retrieval parameters 

        Returns:
            f_grid (numpy array): the frequency grid to do the retrievals.
        """

        # Creation of a cubic decay in the resolution around the observation frequency
        n_f = retrieval_param["number_of_freq_points"]  # Number of points
        bw = 1.3*retrieval_param["bandwidth"]  # Bandwidth
        x = np.linspace(-1, 1, n_f)
        f_grid = x ** 3 + x / retrieval_param["irregularity_f_grid"]
        f_grid = f_grid * bw / (max(f_grid) - min(f_grid)) + \
            retrieval_param['obs_freq']

        #f_grid = np.linspace(retrieval_param["f_min"]-10, retrieval_param["f_max"]+10, n_f)
        #f_grid = np.concatenate((f_grid,np.arange(148.975e9,149.975e9,1e6)))
        #f_grid = np.arange(141e9,150e9,1e6)

        # An option to plot the frequency grid:
        if retrieval_param["show_f_grid"]:
            fig = plt.figure()
            plt.semilogy(f_grid[1:]/1e9, np.diff(f_grid)/1e3, '.')
            plt.xlim((retrieval_param['obs_freq']-200e6) /
                     1e9, (retrieval_param['obs_freq']+200e6)/1e9)
            # plt.ylim(0,300)
            plt.ylabel(r'$\Delta f$ [kHz]')
            plt.suptitle('Frequency grid spacing')
            plt.show()

        return f_grid

    ############################################################################################################################
    # Retrievals grid:
    def set_pgrids(self, retrieval_param):
        """Function to define the pressure grid for ozone and continuum retrievals

        Args:
            retrieval_param (dict): dictionary with all retrieval parameters 

        Returns:
            p_grid_retrieval: the retrieval pressure grid
            p_grid_retrieval_h2o: the retrieval pressure grid for the water vapour continuum
        """ 
        if retrieval_param["retrieval_grid_type"] == 'pressure':
            p_grid_retrieval = np.logspace(5, -1, 61)
        else:
            z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
            z_top_ret = retrieval_param["z_top_ret_grid"]
            z_res_ret = retrieval_param["z_resolution_ret_grid"]
            z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
            #z_grid_retrieval = np.concatenate((np.arange(10e3, 50e3, 2e3),np.arange(50e3, 100e3, 4e3)))
            p_grid_retrieval = z2p_simple(z_grid_retrieval)

        if retrieval_param["retrieval_h2o_grid_type"] == 'pressure':
            if retrieval_param['verbose'] > 1:
                print('Retrieval p_grid for water defined by pressure')
            p_grid_retrieval_h2o = np.array(retrieval_param["h2o_pressure"])
        else:
            if retrieval_param['verbose'] > 1:
                print('Retrieval p_grid from '+str(p_grid_retrieval[0])+' to '+str(p_grid_retrieval[-1]))
            z_bottom_ret_h2o = retrieval_param["z_bottom_ret_grid_h2o"]
            z_top_ret = retrieval_param["z_top_ret_grid_h2o"]
            z_res_ret = retrieval_param["z_resolution_ret_grid_h2o"]
            z_grid_retrieval_h2o = np.arange(z_bottom_ret_h2o, z_top_ret, z_res_ret)
            p_grid_retrieval_h2o = z2p_simple(z_grid_retrieval_h2o)

        return  p_grid_retrieval, p_grid_retrieval_h2o
    
    ############################################################################################################################
    # The retrieval function
    ############################################################################################################################
    def retrieve_cycle(self, spectro_dataset, retrieval_param, ac_sim_FM=None, sensor = None):
        """ Retrieval of a single integration cycle defined in retrieval_param

        Args:
            spectro_dataset (xarray): the daily level 1b dataset
            retrieval_param (dict): dictionary with all retrieval parameters 
            ac_sim_FM (optional): a generic spectra simulated from a previous FM, only for validation purposes. Defaults to None.
            sensor (optional): a sensor object to consider for the observations. Defaults to None. 
                In this case, the object is defined within the retrieval process.

        Raises:
            ValueError: Some error values when retrieval parameters are not correct. TO BE IMPROVED

        Returns:
             ac: an ARTS controler object
             retrieval_param (dict): dictionary with all retrieval parameters 
             sensor: a sensor object
        """
        print('###################################################################################')
        start_time = time.time()
        cycle = retrieval_param["integration_cycle"]
        if ac_sim_FM is None:
            print("Retrieval of ozone and water vapor continuum")
            
            good_channels = spectro_dataset.good_channels[cycle].data == 1
            bad_channels = spectro_dataset.good_channels[cycle].data == 0

            if retrieval_param['use_all_channels']:
                ds_freq = spectro_dataset.frequencies[cycle].values
            else:
                ds_freq = spectro_dataset.frequencies[cycle].values[good_channels]
            
            if retrieval_param['window_corrected_spectrum']:
                print('Using window corrected spectrum interpolated for bad channels')
                ds_y = spectro_dataset.Tb_win_corr[cycle].where(good_channels).interpolate_na(dim='channel_idx', method='nearest', fill_value="extrapolate").values
            else:
                print('Using non window corrected spectrum')
                ds_y = spectro_dataset.Tb[cycle].values[good_channels]
                #ds_y = spectro_dataset.Tb_win_corr[cycle].values[good_channels]
                #ds_y = spectro_dataset.intensity_planck_win_corr[cycle].values[good_channels]  

            ds_num_of_channel = len(ds_freq)
            #ds_Tb = Tb[cycle].values

            ds_bw = max(ds_freq) - min(ds_freq)

            ds_df = ds_bw/(ds_num_of_channel-1)

            retrieval_param["bandwidth"] = self.bandwidth[0]
            # defining simulation grids
            if retrieval_param["binned_ch"]:

                f_grid = self.make_f_grid(retrieval_param)
                ds_y = ds_y[np.arange(0,len(ds_y),18)]
                ds_freq = ds_freq[np.arange(0,len(ds_freq),18)]
                ds_num_of_channel = len(ds_freq)

                ds_df = ds_bw/(ds_num_of_channel-1)

                retrieval_param["bandwidth"] = self.bandwidth[0]
            else:
                if retrieval_param['sensor'] == 'FFT_SB':
                    f_grid = self.make_f_grid_double_sideband(retrieval_param)
                else:
                    f_grid = self.make_f_grid(retrieval_param)

            #print('Minimum of the frequency grid spacing [kHz]: ', min(np.diff(f_grid))/1e3)
        else:
            print("Retrieval of Ozone and H20 providing simulated measurement vector")
            ds_freq = ac_sim_FM.ws.f_backend.value
            ds_y = ac_sim_FM.ws.y.value + np.random.normal(0, 0.2, len(ds_freq)) + 0 + 1e-9*(
                ds_freq-ds_freq[0])*(0)  # Gaussian noise + linear baseline possible
            ds_num_of_channel = len(ds_freq)
            #ds_Tb = Tb[cycle].values

            ds_bw = max(ds_freq) - min(ds_freq)

            ds_df = ds_bw/(ds_num_of_channel-1)

            retrieval_param["bandwidth"] = self.bandwidth[0]+100e6
            f_grid = self.make_f_grid(retrieval_param)

        ############################################################################################ 
        # Iniializing ArtsController object
        ac = arts.ArtsController(verbosity=0, agenda_verbosity=0)
        ac.setup(atmosphere_dim=1, iy_unit='PlanckBT', ppath_lmax=5e3, stokes_dim=1)

        # Storing some geometrical parameters for later: 
        retrieval_param["zenith_angle"] = self.reference_elevation_angle - \
            spectro_dataset.mean_sky_elevation_angle.values[cycle] + retrieval_param['pointing_angle_corr']
        retrieval_param["azimuth_angle"] = spectro_dataset.azimuth_angle.values[cycle]
        retrieval_param["lat"] = spectro_dataset.lat[cycle].values
        retrieval_param["lon"] = spectro_dataset.lon[cycle].values
        retrieval_param["station_altitude"] = spectro_dataset.alt[cycle].values
        
        try:
            retrieval_param["time"] = spectro_dataset.time[cycle].values
            retrieval_param['time_start'] = spectro_dataset.first_sky_time[cycle].values
            retrieval_param['time_stop'] = spectro_dataset.last_sky_time[cycle].values
        except:
            retrieval_param["time"] = 0
            retrieval_param['time_start'] = datetime.date(2019,1,30)
            retrieval_param['time_stop'] = self.date
        
        retrieval_param["f_max"] = max(ds_freq)
        retrieval_param["f_min"] = min(ds_freq)

        print('Time of the measurement: ',pd.to_datetime(retrieval_param["time"]).strftime("%Y-%m-%d, %H:%M:%S"))

        ##################################################
        # Defining simulation grids for the forward model: 
        z_grid = np.arange(retrieval_param["z_bottom_sim_grid"], retrieval_param["z_top_sim_grid"], retrieval_param["z_resolution_sim_grid"])
        p_grid = z2p_simple(z_grid)
        ac.set_grids(f_grid, p_grid)

        ################################################## 
        # Setting spectroscopy
        # abs_species = ["O3","H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252", "O2-PWR93","N2-SelfContStandardType"]
        water_vapor_model = retrieval_param['water_vapor_model']
        ac.set_spectroscopy_from_file(
            abs_lines_file=retrieval_param['line_file'],
            abs_species=retrieval_param['selected_species'],
            format='Arts',
            line_shape=["VVH", 750e9],
        )

        ##################################################
        # Defining atmosphere:
        if retrieval_param['atm'] == 'fascod':
            atm = gromora_atmosphere.get_apriori_fascod(retrieval_param)
            ac.set_atmosphere(atm, vmr_zeropadding=True)
        elif retrieval_param['atm'] == 'ecmwf_cira86':
            ecmwf_prefix = f'ecmwf_oper_v{2}_{self.location}_%Y%m%d.nc'
            retrieval_param['ecmwf_prefix'] = ecmwf_prefix
            try:
                atm = gromora_atmosphere.get_apriori_atmosphere_fascod_ecmwf_cira86(
                    retrieval_param,
                    retrieval_param['ecmwf_store_location'],
                    retrieval_param['cira86_path'],
                    pd.to_datetime(retrieval_param['time_start']),
                    pd.to_datetime(retrieval_param['time_stop']),
                    retrieval_param['extra_time_ecmwf'],
                    z_grid
                )
                ac.set_atmosphere(atm, vmr_zeropadding=True)
            except:
                # Try again another merging method for PTZ profile
                retrieval_param['ptz_merge_method'] = 'max_diff_surf'
                atm = gromora_atmosphere.get_apriori_atmosphere_fascod_ecmwf_cira86(
                    retrieval_param,
                    retrieval_param['ecmwf_store_location'],
                    retrieval_param['cira86_path'],
                    pd.to_datetime(retrieval_param['time_start']),
                    pd.to_datetime(retrieval_param['time_stop']),
                    retrieval_param['extra_time_ecmwf'],
                    z_grid
                )
                ac.set_atmosphere(atm, vmr_zeropadding=True)
        else:
            raise ValueError('Atmosphere type not recognized')
        
        retrieval_param['test_apriori'] = ac.ws.vmr_field.value[0,:,0]

        if retrieval_param["surface_altitude"] < min(ac.ws.z_field.value[:, 0, 0]):
            retrieval_param["surface_altitude"] = min(ac.ws.z_field.value[:, 0, 0])
            retrieval_param["observation_altitude"] = min(ac.ws.z_field.value[:, 0, 0])
            print('Surface altitude has been changed to min of z_field : ',
                  retrieval_param["surface_altitude"])

        ac.set_surface(retrieval_param["surface_altitude"])

        # Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
        # ac.apply_hse(100e2, 0.5)  # value taken from GROMOS retrieval
        ac.apply_hse(ac.ws.p_grid.value[0], 0.5)

        ##################################################
        # Set observation:
        obs = arts.Observation(
            za=retrieval_param["zenith_angle"],
            aa=retrieval_param["azimuth_angle"],
            lat=retrieval_param["lat"],
            lon=retrieval_param["lon"],
            alt=retrieval_param["observation_altitude"],
            time=retrieval_param["time"]
        )

        ac.set_observations([obs])

        ##################################################
        # Define sensor observation:      
        if sensor is None:
            sensor = self.define_sensor(retrieval_param, ds_freq, ds_df)

        ac.set_sensor(sensor)

        ##################################################
        # doing the checks
        ac.checked_calc(negative_vmr_ok=False)

        start_FM_time = time.time()
        if retrieval_param['verbose'] > 0:
            print(f"Setup time: {start_FM_time - start_time:1f} seconds")
        # FM + noise --> to retrieve as test !
        #ac.ws.iy_aux_vars = ["Optical depth"]
        #sprint('Compute Optical Depth')

        if retrieval_param['show_FM']:
            y_FM = ac.y_calc()
            self.plot_FM_comparison(ds_freq, y_FM[0], ds_y)

        if retrieval_param['FM_only']:
            # try to create lookup tables
            # ac.ws.abs_lookupSetupWide(
            #     abs_p=ac.ws.p_grid.value[0],
            # )
           # ac.ws.abs_lookupCalc()
            #ac.ws.jacobianAddBasicCatalogParameter(ac.ws.jacobian_quantities, ac.ws.jacobian_agenda, "O3-666 TR UP J 0 LO J 0", "Line Strength" )

            y_FM = ac.y_calc(jacobian_do=False)
           # plot_FM_comparison(ds_freq, y_FM[0], ds_y)
            return ac, retrieval_param, sensor

        #####################################################################################################################
        # Setup the retrieval
        # Set measurement vector
        ac.set_y([ds_y])

        # if len(ds_y) < 1000:
        #     ac.oem_converged = False
        #     return ac, retrieval_param, sensor
        
        
        ################################################
        # Defining retrieval grids:
        lat_ret_grid = np.array([retrieval_param["lat"]])
        lon_ret_grid = np.array([retrieval_param["lon"]])

        p_grid_retrieval, p_grid_retrieval_h2o = self.set_pgrids(retrieval_param)

        ################################################
        # Setting the apriori covariance matrices in a separate function 
        sx = self.set_o3_apriori_covariance(retrieval_param, p_grid_retrieval)
        sx_water = self.set_continuum_apriori_covariance(retrieval_param, p_grid_retrieval_h2o)

        ################################################
        # Defining retrievals quantities: 
        # Ozone
        ozone_ret = arts.AbsSpecies(
            species='O3',
            p_grid=p_grid_retrieval,
            lat_grid=lat_ret_grid,
            lon_grid=lon_ret_grid,
            covmat=sx,
            unit='vmr',
        )

        # Water vapor continuum
        h2o_ret = arts.AbsSpecies(
            species=water_vapor_model,
            p_grid=p_grid_retrieval_h2o,
            lat_grid=lat_ret_grid,
            lon_grid=lon_ret_grid,
            covmat=sx_water,
            unit='rel',
        )
        ################################################
        # Polyfit
        polyfit_ret = arts.Polyfit(
            poly_order=retrieval_param['poly_order'] , covmats=[np.array([[retrieval_param['covmat_polyfit_0']]]), np.array([[retrieval_param['covmat_polyfit_1']]]), np.array([[retrieval_param['covmat_polyfit_2']]])]
        )
        # polyfit_ret = arts.Polyfit(
        #     poly_order=1, covmats=[np.array([[0.01]]), np.array([[0.1]])]
        # )
        # polyfit_ret = arts.Polyfit(
        #     poly_order=0, covmats=[np.array([[0.01]])]
        # )
        # polyfit_ret = arts.Polyfit(
        #     poly_order=3, covmats=[np.array([[20]]), np.array([[10]]), np.array([[5]]), np.array([[1]])]
        # )

        ################################################
        # Frequency shift
        fshift_ret = arts.FreqShift(100e3, df=50e3)

        ################################################
        # Sinefit

        periods = retrieval_param['sinefit_periods']
        covmat_sinefit = retrieval_param['sinefit_covmat']  
        # covmat_sinefit = covmat.covmat_diagonal_sparse(
        #     np.ones_like([1, 1]))

        # sinefit_ret = arts.RetrievalQuantity(
        #     'Sinefit', covmat_sinefit, period_lengths=periods)
        sinefit_ret = arts.Sinefit(periods, covmats=covmat_sinefit )#  [np.array([1,1]), np.array([1,1])])
        #ac.define_retrieval(retrieval_quantities=[ozone_ret], y_vars=y_var)
        #ac.define_retrieval(retrieval_quantities=[ozone_ret, h2o_ret], y_vars=y_var)

        ##################################################################################################################### 
        # Defining measurement covariance matrix: 
        if ac_sim_FM is None:
            if retrieval_param['noise_covariance']  == 'noise_level':
                y_var = retrieval_param['increased_var_factor']*np.square(
                    spectro_dataset.noise_level[cycle].data) * np.ones_like(ds_y)
                #y_var = retrieval_param['increased_var_factor']*np.square(np.std(np.diff(ds_y)/np.sqrt(2))) * np.ones_like(ds_y)
                if len(y_var) == len(bad_channels):
                    # increase variance for spurious channels by some factor
                    y_var[bad_channels] = 1e5*np.square(spectro_dataset.noise_level[cycle].data)

            elif retrieval_param['noise_covariance']  == 'stdTb':
                y_var = retrieval_param['increased_var_factor']*np.square(
                    spectro_dataset.mean_std_Tb[cycle].data) * np.ones_like(ds_y)
            else:
                raise ValueError('Please select an existing noise covariance type !')
        else:
            # For cases with simulated measurement vector 
            print('Using standard y var')
            y_var = 0.04 * np.ones_like(ds_y)

        if retrieval_param['verbose'] > 0:
            print(f'Noise level for meas. cov: {np.sqrt(np.median(y_var)):.2f} K')
        
        ac.noise_variance_vector = y_var
        ac.tropospheric_opacity = spectro_dataset.tropospheric_opacity[cycle].values

        ##################################################################################################################### 
        # Apply the retrievals quantites defined by the user: 
        if retrieval_param['retrieval_quantities'] == 'o3_h2o':
            print('Retrievals quantities: O3 and H2O continuum')
            ac.define_retrieval(retrieval_quantities=[
                                ozone_ret, h2o_ret], y_vars=y_var)
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift':
            print('Retrievals quantities: O3, H2O continuum and Fshift')
            ac.define_retrieval(retrieval_quantities=[
                                ozone_ret, h2o_ret, fshift_ret],  y_vars=y_var)
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_polyfit':
            print('Retrievals quantities: O3, H2O continuum and Polyfit')
            ac.define_retrieval(retrieval_quantities=[
                                ozone_ret, h2o_ret, polyfit_ret],  y_vars=y_var)
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit':
            print('Retrievals quantities: O3, H2O continuum, Fshift and Polyfit')
            ac.define_retrieval(retrieval_quantities=[
                                ozone_ret, h2o_ret, polyfit_ret, fshift_ret],  y_vars=y_var)
        elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit_sinefit':
            print('Retrievals quantities: O3, H2O continuum, Fshift, Polyfit and Sinefit')
            ac.define_retrieval(retrieval_quantities=[
                                ozone_ret, h2o_ret, polyfit_ret, fshift_ret, sinefit_ret],  y_vars=y_var)
        else:
            raise ValueError('Retrieval quantities not recognized !')

        retrievals_setup_time = time.time()
        if retrieval_param['verbose'] > 1:
            print(f"retrievals_setup_time: {retrievals_setup_time - start_FM_time:.1} s")

        # Run retrieval (parameter taken from MOPI)
        # SOMORA is using 'lm': Levenberg-Marquardt (LM) method

        #####################################################################################################################
        # Agenda for GROMORA retrievals
        # see: https://atmtools.github.io/arts-docs-master/docserver/agendas/inversion_iterate_agenda.html
        @arts_agenda
        def gromora_inversion_agenda(ws):
            """Custom inversion iterate agenda to ignore bad partition functions."""
            ws.Ignore(ws.inversion_iteration_counter)

            ws.xClip(ijq=0, limit_low=0.00000000001, limit_high=0.00002)

            # Map x to ARTS' variables
            ws.x2artsAtmAndSurf()
            ws.x2artsSensor()

            # To be safe, rerun some checkss
            ws.atmfields_checkedCalc(negative_vmr_ok=True)
            ws.atmgeom_checkedCalc()

            # Calculate yf and Jacobian matching x
            ws.yCalc() #()ws.yf

            # Add baseline term
            #ws.VectorAddElementwise(ws.yf, ws.y, ws.y_baseline)
            ws.VectorAddVector(ws.yf, ws.y, ws.y_baseline)

            # This method takes cares of some "fixes" that are needed to get the Jacobian
            # right for iterative solutions. No need to call this WSM for linear inversions.
            ws.jacobianAdjustAndTransform()

        #####################################################################################################################
        # Run the optimal estimation method in ARTS:
        ac.checked_calc()
        ac.oem(
            method=retrieval_param['oem_method'],
            max_iter=retrieval_param['max_iter'],
            stop_dx=retrieval_param['stop_dx'],
            display_progress=1,
            lm_ga_settings=retrieval_param['lm_ga_setting'] ,
            inversion_iterate_agenda=gromora_inversion_agenda,
        )

        retrievals_time = time.time()
        if retrieval_param['verbose'] > 0:
            print(f"Total retrieval time: {retrievals_time - retrievals_setup_time:.2f} seconds")

            print(f"Convergence status ORM (0 = converged): {ac.oem_diagnostics[0]:.1f}")
            print(f"End value of the cost function: {ac.oem_diagnostics[2]:.4f}")
            print(f"Number of iteration used: {ac.oem_diagnostics[4]:.1f}")

        print('#####################################')

        if not ac.oem_converged:
            print("OEM did not converge.")
            print("OEM diagnostics: " + str(ac.oem_diagnostics))
            for e in ac.oem_errors:
                print("OEM error: " + e)
                continue
    
        return ac, retrieval_param, sensor

    #####################################################################################################################
    # Level 2 writing function
    #####################################################################################################################
    def write_level2_gromora(self, level2, retrieval_param, full_name):
        """Function to format and write level2 from GROMOS and SOMORA into netCDF files.

        Args:
            level2 (xarray dataset): the retrieved level 2
            retrieval_param (dict): dictionary with all retrieval parameters 
            full_name (str): the full name of the output file
        """
        ##################################################
        # Delete some useless variable
        level2 = level2.isel(o3_lat=0,o3_lon=0, h2o_continuum_lat=0, h2o_continuum_lon=0, drop=True)# .drop_vars(['o3_lat','o3_lon', 'h2o_continuum_lat', 'h2o_continuum_lon']).reset_coords(names=['o3_lat','o3_lon', 'h2o_continuum_lat', 'h2o_continuum_lon'], drop=True )
        level2 = level2.drop_vars(['obs_time'])

        ##################################################
        # Time coordinates_ 
        level2.time.attrs['standard_name'] = 'time'
        # level2.time.attrs['calendar'] = 'proleptic_gregorian'
        # level2.time.attrs['units'] = 'days since 2000-01-01 00:00:00'
        level2.time.encoding['units'] = 'days since 2000-01-01 00:00:00'
        level2.time.encoding['calendar'] = 'proleptic_gregorian'
        level2.time.attrs['timezone'] = 'Z'
        level2.time.attrs['description'] = 'mean time recorded at the beginning of all sky measurements during this integration cycle'

        # adding local solar time and MJD2K 
        julian_dates = mjd2k_date(pd.to_datetime(level2.time.data))

        local_solar_time = list()
        solar_zenith_angle = list()
        for t in level2.time.values:
            lst, ha, sza, night, tc = get_LST_from_GROMORA(t, retrieval_param['lat'], retrieval_param['lon'])
            local_solar_time.append(lst)
            solar_zenith_angle.append(sza)
        
        level2['local_solar_time'] =  ('time', local_solar_time)
        level2.local_solar_time.encoding['calendar'] = 'proleptic_gregorian'
        level2.local_solar_time.encoding['units'] = 'days since 2000-01-01 00:00:00'
        level2.local_solar_time.attrs['description'] = 'local solar time computed from the mean measurement time'
        
        level2['MJD2K'] =  ('time', julian_dates)
        level2.MJD2K.attrs['standard_name'] = 'MJD2K'
        level2.MJD2K.attrs['long_name'] = 'Modified Julian Date 2000'
        level2.MJD2K.attrs['units'] = 'MJD2K'
        level2.MJD2K.attrs['description'] = 'MJD2K as defined by GEOMS: it is 0.000000 on January 1, 2000 at 00:00:00 UTC'

        ##################################################
        # Global attributes 
        level2.attrs['title'] = 'Ozone retrievals from microwave ground-based radiometer'
        level2.attrs['location'] = self.location
        level2.attrs['instrument'] = self.instrument_name
        level2.attrs['source'] = self.source
        level2.attrs['name'] = self.name_PI
        level2.attrs['featureType'] = 'timeSeries'
        level2.attrs['contact'] = self.contact
        level2.attrs['references'] = ''
        level2.attrs['comment'] = ''
        level2.attrs['institution'] = self.institution
        level2.attrs['number_of_spectrometer'] = self.number_of_spectrometer 
        level2.attrs['filename_level1a'] = self.filename_level1a
        level2.attrs['filename_level1b'] = self.filename_level1b
        level2.attrs['raw_data_software_version'] = self.raw_data_software_version
        level2.attrs['calibration_version'] = self.calibration_version

        level2.attrs['spectroscopy'] = retrieval_param['line_file']

        level2.attrs['outlier_detection'] = self.outlier_detection

        level2.attrs['raw_data_filename'] = self.raw_data_filename
        level2.attrs['raw_data_software_version'] =  self.raw_data_software_version
        # level2.attrs['raw_file_warning'] = self.raw_file_warning
        level2.attrs['labview_logfile_warning'] = self.labview_logfile_warning
        level2.attrs['filtering_of_calibrated_spectra'] = self.filtering_of_calibrated_spectra

        level2.attrs['data_start_date'] = np.datetime_as_string(level2.time.data[0], timezone='UTC')
        level2.attrs['data_stop_date'] = np.datetime_as_string(level2.time.data[-1], timezone='UTC')

        level2.attrs['integration_time'] = self.int_time

        # Some extra infos about the retrievals 
        level2.attrs['atmosphere'] = retrieval_param['atm']
        level2.attrs['sensor'] = retrieval_param['sensor'] 
        level2.attrs['retrieved_quantities'] = retrieval_param['retrieval_quantities']
        level2.attrs['spectroscopy'] = retrieval_param['line_file']

        ##################################################
        # Adding atttributes to coordinates:
        level2.o3_p.attrs['standard_name'] = 'o3_p'
        level2.o3_p.attrs['long_name'] = 'retrieval pressure grid for ozone'
        level2.o3_p.attrs['units'] = 'Pa'
        level2.o3_p.attrs['description'] = 'Pressure grid for ozone retrievals'

        level2.h2o_continuum_p.attrs['standard_name'] = 'h2o_continuum_p'
        level2.h2o_continuum_p.attrs['long_name'] = 'retrieval pressure grid for continuum'
        level2.h2o_continuum_p.attrs['units'] = 'Pa'
        level2.h2o_continuum_p.attrs['description'] = 'Pressure grid for continuum retrievals'

        level2.f.attrs['standard_name'] = 'frequency'
        level2.f.attrs['long_name'] = 'measurement_frequency'
        level2.f.attrs['units'] = 'Hz'
        level2.f.attrs['description'] = 'measurement frequency grid'

        ##################################################
        # Adding atttributes to important variables:
        # Ozone
        level2.o3_x.attrs['standard_name'] = 'ozone_profile'
        level2.o3_x.attrs['long_name'] = 'ozone profile volume mixing ratio'
        level2.o3_x.attrs['units'] = 'VMR'
        level2.o3_x.attrs['description'] = 'Ozone profile retrieved with OEM'
        level2.o3_x.attrs['valid_min'] = 0
        level2.o3_x.attrs['valid_max'] = 50e-6

        level2.o3_xa.attrs['standard_name'] = 'o3_xa'
        level2.o3_xa.attrs['long_name'] = 'ozone apriori'
        level2.o3_xa.attrs['units'] = 'VMR'
        level2.o3_xa.attrs['description'] = retrieval_param['o3_apriori']
        level2.o3_xa.attrs['covariance'] = retrieval_param['o3_apriori_covariance']

        level2.o3_mr.attrs['standard_name'] = 'o3_mr'
        level2.o3_mr.attrs['long_name'] = 'ozone measurement response'
        level2.o3_mr.attrs['units'] = '1'

        level2.o3_eo.attrs['standard_name'] = 'o3_eo'
        level2.o3_eo.attrs['long_name'] = 'observation error for ozone'
        level2.o3_eo.attrs['units'] = 'VMR'

        level2.o3_es.attrs['standard_name'] = 'o3_es'
        level2.o3_es.attrs['long_name'] = 'smoothing error for ozone'
        level2.o3_es.attrs['units'] = 'VMR'
        
        level2.o3_avkm.attrs['standard_name'] = 'o3_avkm'
        level2.o3_avkm.attrs['long_name'] = 'averaging kernels for ozone'
        level2.o3_avkm.attrs['units'] = '1'

        level2.o3_z.attrs['standard_name'] = 'o3_z'
        level2.o3_z.attrs['long_name'] = 'geometric altitude for ozone'
        level2.o3_z.attrs['units'] = 'm'
    
        level2.o3_fwhm.attrs['standard_name'] = 'o3_fwhm'
        level2.o3_fwhm.attrs['long_name'] = 'full width at half maximum for ozone'
        level2.o3_fwhm.attrs['units'] = 'm'

        level2.o3_offset.attrs['standard_name'] = 'o3_offset'
        level2.o3_offset.attrs['long_name'] = 'altitude offset for ozone AVKs'
        level2.o3_offset.attrs['units'] = 'm'
    
        # Continuum
        level2.h2o_continuum_x.attrs['standard_name'] = 'continuum_retrieval'
        level2.h2o_continuum_x.attrs['long_name'] = 'water vapor continuum'
        level2.h2o_continuum_x.attrs['units'] = '1'
        level2.h2o_continuum_x.attrs['continuum model'] =  retrieval_param['water_vapor_model']
        level2.h2o_continuum_x.attrs['description'] = 'ratio to apriori value'

        level2.h2o_continuum_xa.attrs['standard_name'] = 'h2o_continuum_xa'
        level2.h2o_continuum_xa.attrs['long_name'] = 'water vapor continuum apriori'
        level2.h2o_continuum_xa.attrs['units'] = '1'

        level2.h2o_continuum_mr.attrs['standard_name'] = 'h2o_continuum_mr'
        level2.h2o_continuum_mr.attrs['long_name'] = 'continuum measurement response'
        level2.h2o_continuum_mr.attrs['units'] = '1'

        level2.h2o_continuum_eo.attrs['standard_name'] = 'h2o_continuum_eo'
        level2.h2o_continuum_eo.attrs['long_name'] = 'observation error for continuum'
        level2.h2o_continuum_eo.attrs['units'] = 'VMR'

        level2.h2o_continuum_es.attrs['standard_name'] = 'h2o_continuum_es'
        level2.h2o_continuum_es.attrs['long_name'] = 'smoothing error for continuum'
        level2.h2o_continuum_es.attrs['units'] = 'VMR'
        
        level2.h2o_continuum_avkm.attrs['standard_name'] = 'h2o_continuum_avkm'
        level2.h2o_continuum_avkm.attrs['long_name'] = 'averaging kernels for continuum'
        level2.h2o_continuum_avkm.attrs['units'] = '1'

        if 'temperature_profile' in list(level2.keys()):
            level2.temperature_profile.attrs['standard_name'] = 'temperature_profile'
            level2.temperature_profile.attrs['long_name'] = 'atmospheric temperature profiles interpolated on retrieval grid'
            level2.temperature_profile.attrs['units'] = 'K'

        # Polyfit
        if 'poly_fit_x' in list(level2.keys()):
            
            level2.poly_order.attrs['standard_name'] = 'poly_order'
            level2.poly_order.attrs['long_name'] = 'poly_order'
            level2.poly_order.attrs['units'] = '1'
            level2.poly_order.attrs['description'] = 'degree order for polyfit'

            level2.poly_fit_x.attrs['standard_name'] = 'polyfit_retrieval'
            level2.poly_fit_x.attrs['long_name'] = 'polyfit retrievals'
            level2.poly_fit_x.attrs['units'] = 'K'
            level2.poly_fit_x.attrs['poly_order'] = retrieval_param['poly_order']
            level2.poly_fit_x.attrs['description'] = ''

            level2.poly_fit_xa.attrs['standard_name'] = 'poly_fit_xa'
            level2.poly_fit_xa.attrs['long_name'] = 'polyfit apriori'
            level2.poly_fit_xa.attrs['units'] = '1'

            level2.poly_fit_mr.attrs['standard_name'] = 'poly_fit_mr'
            level2.poly_fit_mr.attrs['long_name'] = 'polyfit measurement response'
            level2.poly_fit_mr.attrs['units'] = '1'

        # Frequency shift
        if 'freq_shift_x' in list(level2.keys()):
            level2.f_shift_grid.attrs['standard_name'] = 'f_shift_grid'
            level2.f_shift_grid.attrs['long_name'] = 'f_shift_grid'
            level2.f_shift_grid.attrs['units'] = '1'
            level2.f_shift_grid.attrs['description'] = ''

            level2.freq_shift_x.attrs['standard_name'] = 'freq_shift_retrieval'
            level2.freq_shift_x.attrs['long_name'] = 'frequency shift retrievals'
            level2.freq_shift_x.attrs['units'] = 'Hz'
            level2.freq_shift_x.attrs['description'] = ''

            level2.freq_shift_xa.attrs['standard_name'] = 'freq_shift_xa'
            level2.freq_shift_xa.attrs['long_name'] = 'frequency shift apriori'
            level2.freq_shift_xa.attrs['units'] = 'Hz'

            level2.freq_shift_mr.attrs['standard_name'] = 'freq_shift_mr'
            level2.freq_shift_mr.attrs['long_name'] = 'frequency shift measurement response'
            level2.freq_shift_mr.attrs['units'] = '1'

        level2.y.attrs['standard_name'] = 'brightness_temperature'
        level2.y.attrs['long_name'] = 'measurement vector'
        level2.y.attrs['units'] = 'K'
        level2.y.attrs['description'] = 'integrated brightness temperature for this cycle, corrected for window and troposphere'

        level2.bad_channels.attrs['standard_name'] = 'bad_channels'
        level2.bad_channels.attrs['long_name'] = 'bad channels identified on measurement vector'
        level2.bad_channels.attrs['units'] = '-'
        level2.bad_channels.attrs['description'] = 'a boolean vector identifying the bad channels on the measurement vector'

        level2.yf.attrs['standard_name'] = 'fitted_measurement'
        level2.yf.attrs['long_name'] = 'fitted measurement vector'
        level2.yf.attrs['units'] = 'K'
        level2.yf.attrs['description'] = 'fitted window corrected brightness temperature'

        level2.median_noise.attrs['standard_name'] = 'noise_level'
        level2.median_noise.attrs['long_name'] = 'radiometric noise level'
        level2.median_noise.attrs['units'] = 'K'
        level2.median_noise.attrs['description'] = 'noise level used as measurement error'

        level2.y_baseline.attrs['standard_name'] = 'measurement_baseline'
        level2.y_baseline.attrs['long_name'] = 'baseline measurement vector '
        level2.y_baseline.attrs['units'] = 'K'
        level2.y_baseline.attrs['description'] = 'retrieved baseline for measurement vector'
        
        level2.obs_za.attrs['standard_name'] = 'sensor_zentih_angle'
        level2.obs_za.attrs['long_name'] = 'zenith angle'
        level2.obs_za.attrs['units'] = 'degree'

        level2.obs_aa.attrs['standard_name'] = 'sensor_azimuth_angle'
        level2.obs_aa.attrs['long_name'] = 'azimuth angle'
        level2.obs_aa.attrs['units'] = 'degree'
        level2.obs_aa.attrs['description'] = 'angle measured clockwise positive, 0 deg is northwise'

        if 'tropospheric_opacity' in list(level2.keys()):
            level2.tropospheric_opacity.attrs['standard_name'] = 'tropospheric_opacity'
            level2.tropospheric_opacity.attrs['long_name'] = 'tropospheric_opacity computed with Ingold method during calibration'
            level2.tropospheric_opacity.attrs['units'] = 'Np'

        level2['solar_zenith_angle'] =  ('time', solar_zenith_angle)
        level2.solar_zenith_angle.attrs['standard_name'] = 'solar_zenith_angle'
        level2.solar_zenith_angle.attrs['long_name'] = 'solar zenith angle'
        level2.solar_zenith_angle.attrs['units'] = 'deg'
        level2.solar_zenith_angle.attrs['description'] = 'angle between the sun rays and zenith, minimal at local solar noon'

        ##################################################
        # Renaming some variable:
        level2 = level2.rename_vars( {'obs_lat':'lat','obs_lon':'lon', 'obs_alt':'alt', } ) 
        
        level2.lat.attrs['standard_name'] = 'latitude'
        level2.lat.attrs['long_name'] = 'station latitude'
        level2.lat.attrs['units'] = 'degree_north'
        level2.lat.attrs['description'] = 'latitude defined according to WGS84'

        level2.lon.attrs['standard_name'] = 'longitude'
        level2.lon.attrs['long_name'] = 'station longitude'
        level2.lon.attrs['units'] = 'degree_east'
        level2.lon.attrs['description'] = 'longitude defined according to WGS84'

        level2.alt.attrs['standard_name'] = 'altitude'
        level2.alt.attrs['long_name'] = 'station altitude'
        level2.alt.attrs['units'] = 'm'
        level2.alt.attrs['description'] = 'altitude above see level'

        level2.oem_diagnostics.attrs['standard_name'] = 'oem_diagnostics'
        level2.oem_diagnostics.attrs['long_name'] = 'optimal estimation diagnostics'
        level2.oem_diagnostics.attrs['units'] = '1'
        level2.oem_diagnostics.attrs['description'] = 'Vector of retrievals outputs from ARTS used for diagnostics'
        level2.oem_diagnostics.attrs['diagnosticValue0'] = 'Convergence status, with coding (0 = converged).'
        level2.oem_diagnostics.attrs['diagnosticValue1'] = 'Start value of cost function.'
        level2.oem_diagnostics.attrs['diagnosticValue2'] = 'End value of cost function.'
        level2.oem_diagnostics.attrs['diagnosticValue3'] = 'End value of y-part of cost function.'
        level2.oem_diagnostics.attrs['diagnosticValue4'] = 'Number of iterations used.'

        # Polyfit
        if 'sine_grid' in list(level2.coords.keys()):
            level2.sine_grid.attrs['standard_name'] = 'sine_grid'
            level2.sine_grid.attrs['long_name'] = 'sine_grid'
            level2.sine_grid.attrs['units'] = '1'
            level2.sine_grid.attrs['description'] = 'the grid for the sinefit retrievals, first element is sine and second is cosine term'

            for i, p in enumerate(retrieval_param['sinefit_periods']):
                p_MHz = f'{p/1e6:.0f}'
                varname = f'sine_fit_{i}_x'
                level2[varname].attrs['standard_name'] = varname
                level2[varname].attrs['long_name'] = 'sinefit_retrieval'
                level2[varname].attrs['units'] = '1'
                level2[varname].attrs['period MHz'] = p/1e6
                level2[varname].attrs['description'] = 'Sinusoidal baseline retrieved for this period, first element is sine and second is cosine term'

        ##################################################
        # Saving it as netCDF v4
        level2.to_netcdf(path=full_name, format='NETCDF4', unlimited_dims='time')

        return level2


    #####################################################################################################################
    # Plotting diagnostics for level 2
    #####################################################################################################################
    def plot_ozone_sel(self, level2_data, outName, spectro, cycles=None, altitude = False, add_baselines=False, to_ppm=1):
        """Plotting function for the retrieved ozone profile

        Args:
            level2_data (xarray dataset): level 2 dataset with all diagnostic quantities
            outName (str): name for the output file
            spectro (str): the spectrometer to use (mostly for MOPI)
            cycles (list, optional): list with cycle number to plot. Defaults to None -> plot all cycle in level 2
            altitude (bool, optional): add altitude coordinate. Defaults to False.
            add_baselines (bool, optional): plot the baseline on the residual plot. Defaults to False.
            to_ppm (int, optional): conversion factor to plot in PPMV. Defaults to 1.
        """
        F0 = self.observation_frequency

        col = self.basecolor

        figure_o3_sel=list()
        fs=28

        if cycles is None:
            cycles = np.arange(len(level2_data[spectro].time))

        for i in cycles:
            f_backend = np.where(level2_data[spectro].bad_channels[i]==0,level2_data[spectro].f.data,  np.nan)
            y = np.where(level2_data[spectro].bad_channels[i]==0, level2_data[spectro].y[i].data, np.nan)# level2_data[spectro].y[i].data
            yf = np.where(level2_data[spectro].bad_channels[i]==0, level2_data[spectro].yf[i].data, np.nan)# level2_data[spectro].yf[i].data
            bl = np.where(level2_data[spectro].bad_channels[i]==0, level2_data[spectro].y_baseline[i].data, np.nan)# level2_data[spectro].y_baseline[i].data 
            r = y - yf
            r = np.where(np.isnan(r), 0, r)
            fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(18,12))
            if self.instrument_name == 'GROMOS':
                print('Binning spectra for similar resolution as SOMORA')
                r = np.convolve(r, np.ones(2) / 2, mode="same")
                axs[0].plot((np.convolve(f_backend, np.ones(2)/2, mode='same') - F0) / 1e6, np.convolve(y, np.ones(2)/2, mode='same'), color='silver', label="observed", alpha=.75)
            else:
                axs[0].plot((f_backend - F0) / 1e6, y, color='silver', label="observed", alpha=.75)
            r_smooth = np.convolve(r, np.ones(128) / 128, mode="same")
            axs[0].plot((f_backend - F0) / 1e6, yf, color='k', label="fitted")
            axs[0].set_ylabel("$T_B$ [K]", fontsize=fs)
            axs[0].set_ylim(np.nanmedian(yf)-4, np.nanmedian(yf)+20)
           # axs[0].set_xlim(-0.5,15)
            axs[0].legend(fontsize=fs)
            axs[1].plot((f_backend - F0) / 1e6, r, color='silver', label="residuals", alpha=.75)
            axs[1].plot((f_backend - F0) / 1e6, r_smooth, color='k', label="residuals smooth")
            if add_baselines:
                axs[1].plot((f_backend - F0) / 1e6, bl, label="baseline")

            axs[1].set_ylim(-4, 4)
            axs[1].legend(fontsize=fs)
            axs[1].set_xlabel("RF - {:.3f} GHz [MHz]".format(F0 / 1e9), fontsize=fs)
            axs[1].set_ylabel(r"$\Delta T_B$ [K]", fontsize=fs)
            for ax in axs:
                ax.grid()
                ax.set_xlim([np.nanmin((f_backend - F0) / 1e6), np.nanmax((f_backend - F0) / 1e6)])
                ax.tick_params(axis='both', which='major', labelsize=fs)
            fig.suptitle(self.instrument_name+' O$_3$ spectrum: '+pd.to_datetime(level2_data[spectro].time[i].data).strftime('%Y-%m-%d %H:%M'), fontsize=fs+4)

            fig.tight_layout(rect=[0, 0.03, 1, 0.99])
            figure_o3_sel.append(fig)


            fig, axs = plt.subplots(nrows=1, ncols=4, sharey=True, figsize=(24,16))

            o3 = level2_data[spectro].isel(time=i).o3_x
            o3_apriori = level2_data[spectro].isel(time=i).o3_xa
            o3_z = level2_data[spectro].isel(time=i).o3_z
            fwhm=level2_data[spectro].isel(time=i).o3_fwhm 
            offset=level2_data[spectro].isel(time=i).o3_offset
            o3_p = level2_data[spectro].isel(time=i).o3_p
            mr = level2_data[spectro].isel(time=i).o3_mr
            #error = lvl2[spectro].isel(time=i).o3_eo +  lvl2[spectro].isel(time=i).o3_es
            error = np.sqrt(level2_data[spectro].isel(time=i).o3_eo**2 +  level2_data[spectro].isel(time=i).o3_es**2)
            error_frac = error/o3
            o3_good = o3.where(mr>0.8).data

            if altitude:
                y_axis=o3_z/1e3
                y_lab = 'Altitude [km]'
            else:
                y_axis = o3_p/100
                y_lab = 'Pressure [hPa] '
            axs[0].fill_betweenx(y_axis, (o3-error)*to_ppm,(o3+error)*to_ppm, color=col, alpha=0.5)
            axs[0].plot(o3*to_ppm, y_axis,'-', linewidth=1.5, label='retrieved',color=col)
            axs[0].plot(o3_apriori*to_ppm, y_axis, '--', linewidth=1.5, label='apriori',color='k')
            #axs[0].set_title('O$_3$ VMR')
            axs[0].set_xlim(-0.5,9)
            if altitude:
                axs[0].set_ylim(5,85)
                axs[0].yaxis.set_major_locator(MultipleLocator(10))
                axs[0].yaxis.set_minor_locator(MultipleLocator(5))
            else:
                axs[0].set_yscale('log')
                axs[0].invert_yaxis()
                axs[0].set_ylim(500,0.005)
               # axs[0].yaxis.set_major_locator(MultipleLocator(10))
              #  axs[0].yaxis.set_minor_locator(MultipleLocator(5))
                axs[0].yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))
            axs[0].set_xlabel('O$_3$ VMR [ppmv]', fontsize=fs)

            axs[0].xaxis.set_major_locator(MultipleLocator(4))
            axs[0].xaxis.set_minor_locator(MultipleLocator(1))
            axs[0].grid(which='both',  axis='x', linewidth=0.5)
            axs[0].set_ylabel(y_lab, fontsize=fs)
            axs[0].legend(fontsize=fs)
            axs[1].plot(mr/4, y_axis,color='k', label='MR/4')
            axs[0].text(
            0.9,
            0.01,
            'a)',
            transform=axs[0].transAxes,
            verticalalignment="bottom",
            horizontalalignment="left",
            fontsize=fs+8
            )
            counter=0
            color_count = 0
            for j, avk in enumerate(level2_data[spectro].isel(time=i).o3_avkm):
                if 0.6 <= np.sum(avk) <= 1.4:
                    counter=counter+1
                    if np.mod(counter,8)==0:
                        axs[1].plot(avk, y_axis, color=cmap(color_count*0.25+0.01))#label='z = '+f'{o3_z.sel(o3_p=avk.o3_p).values/1e3:.0f}'+' km'
                        color_count = color_count +1
                    else:
                        if counter==1:
                            axs[1].plot(avk,y_axis, color='silver', label='AVKs')
                        else:
                            axs[1].plot(avk, y_axis, color='silver')


            # counter=0
            # for avk in level2_data[spectro].isel(time=i).o3_avkm:
            #     if 0.8 <= np.sum(avk) <= 1.2:
            #         counter=counter+1
            #         if np.mod(counter,5)==0:
            #             axs[1].plot(avk, o3_z / 1e3, label='z='+f'{o3_z.sel(o3_p=avk.o3_p).values/1e3:.0f}'+'km', color='r')
            #         else:
            #             axs[1].plot(avk, o3_z / 1e3, color='k')
            axs[1].set_xlabel("Averaging Kernels", fontsize=fs)
            axs[1].set_ylabel("", fontsize=fs)
            axs[1].set_xlim(-0.08,0.4)
            axs[1].xaxis.set_major_locator(MultipleLocator(0.2))
            axs[1].xaxis.set_minor_locator(MultipleLocator(0.05))
            axs[1].legend(loc=1, fontsize=fs-2)
            axs[1].grid(which='both',  axis='x', linewidth=0.5)
            axs[1].text(
            0.9,
            0.01,
            'b)',
            transform=axs[1].transAxes,
            verticalalignment="bottom",
            horizontalalignment="left",
            fontsize=fs+8
            )


            axs[2].plot(level2_data[spectro].isel(time=i).o3_es * 1e6, y_axis, '-', color='k', label="smoothing error")
            axs[2].plot(level2_data[spectro].isel(time=i).o3_eo * 1e6, y_axis, '--' ,color='k', label="measurement error")
            axs[2].set_xlabel("Errors [ppmv]", fontsize=fs)
            axs[2].set_ylabel("", fontsize=fs)
            axs[2].set_xlim(-0.08,1)
            axs[2].xaxis.set_major_locator(MultipleLocator(0.5))
            axs[2].xaxis.set_minor_locator(MultipleLocator(0.1))
            axs[2].legend(loc=1, fontsize=fs-2)
            axs[2].grid(axis='x', linewidth=0.5)
            axs[2].text(
            0.9,
            0.01,
            'c)',
            transform=axs[2].transAxes,
            verticalalignment="bottom",
            horizontalalignment="left",
            fontsize=fs+8
        )

            axs[3].plot(fwhm/1e3, y_axis, color='k', label='FWHM')
            axs[3].plot(offset/1e3, y_axis, '--', color='k', label='AVKs offset')
            axs[3].set_xlim(-15,20)
            axs[3].set_xlabel("Resolution and offset [km]", fontsize=fs)
            axs[3].set_ylabel("", fontsize=fs)
            axs[3].xaxis.set_minor_locator(MultipleLocator(5))
            axs[3].xaxis.set_major_locator(MultipleLocator(10))
            axs[3].grid(which='both', axis='x', linewidth=0.5)
            axs[3].legend(loc=2, fontsize=fs-2)
            axs[3].text(
            0.9,
            0.01,
            'd)',
            transform=axs[3].transAxes,
            verticalalignment="bottom",
            horizontalalignment="left",
            fontsize=fs+8
        )

            #axs[3].plot(level2_data[spectro].isel(time=i).h2o_x * 1e6, o3_z / 1e3, label="retrieved")
            # axs[3].set_xlabel("$VMR$ [ppm]")
            # axs[3].set_ylabel("Altitude [km]")
            # axs[3].legend()
            #axs[3].grid(axis='x', linewidth=0.5)

            # adding altitude axis, thanks Leonie :)
            y1z=1e-3*o3_z.sel(o3_p=48696, tolerance=100,method='nearest')
            y2z=1e-3*o3_z.sel(o3_p=0.5, tolerance=1,method='nearest')
            ax2 = axs[3].twinx()
            ax2.set_yticks(level2_data[spectro].isel(time=i).o3_z) #ax2.set_yticks(altitude)
            ax2.set_ylim(y1z,y2z)
            fmt = FormatStrFormatter("%.0f")
            loc=MultipleLocator(base=10)
            ax2.yaxis.set_major_formatter(fmt)
            ax2.yaxis.set_major_locator(loc)
            ax2.set_ylabel('Altitude [km] ', fontsize=fs)
            ax2.tick_params(axis='both', which='major', labelsize=fs)

            for a in axs:
                #a.set_ylim(10,80)
                a.grid(which='both', axis='y', linewidth=0.5)
                a.grid(which='both', axis='x', linewidth=0.5)
                a.tick_params(axis='both', which='major', labelsize=fs)
            fig.suptitle('O$_3$ retrievals: '+pd.to_datetime(level2_data[spectro].time[i].data).strftime('%Y-%m-%d %H:%M'), fontsize=fs+4)
            fig.tight_layout(rect=[0, 0.01, 1, 0.99])
            figure_o3_sel.append(fig)

        save_single_pdf(outName+'.pdf',figure_o3_sel)    

    def plot_level1b_TB(self, title='', save=False, outfolder='', save_name='int_spectra', idx=None):
        figures = list()
        
        if idx is None:
            figures.append(GROMORA_library.plot_Tb_all(self, self.integrated_data, title=title)) 
        else:
            figures.append(GROMORA_library.plot_Tb_selected(self, self.integrated_data, title=title, idx=idx)) 

        if save:
            save_single_pdf(outfolder+self.instrument_name+'/'+save_name+self.datestr+'.pdf', figures)
            print('saved in '+outfolder+self.instrument_name+'/'+save_name+self.datestr+'.pdf')
            #save_pngs(self.level1_folder+'time_series_'+self.datestr+'_', figures)

    def plot_FM_comparison(self, ds_freq, y_FM, y_obs):
        '''
        short function plotting the result of the FM vs the observation
        '''
        fig = plt.figure()

        fig.suptitle('Comparison between FM and Observation')
        ax = fig.add_subplot(111)
        ax.plot(ds_freq/1e9, y_obs, 'r')
        ax.plot(ds_freq/1e9, y_FM, 'b-', linewidth=2)
        ax.set_ylim((10,250))

        ax.set_xlabel('freq [GHz]')
        ax.set_ylabel('Tb [K]')

        ax.legend(('Observed', 'FM'))

        plt.show()
        pass

    def retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, f_bin = None, tb_bin = None, ac=None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False
 
        retrieval_param['ref_elevation_angle'] = 90

        return retrieval_module.retrieve_cycle_tropospheric_corrected(self, spectro_dataset, retrieval_param, ac_sim_FM=ac)

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
        return GROMORA_library.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    
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
        return GROMORA_library.smooth_corr_spectra(level1b_dataset, retrieval_param)
    
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
        return GROMORA_library.find_bad_channels(self.level1b_ds, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
    
    def find_bad_channels_stdTb(self, spectrometers, stdTb_threshold, apply_on='cal', dimension=['time','channel_idx']):
        """Identification of bad channels for GROMORA retrievals

        Args:
            spectrometers (_type_): _description_
            stdTb_threshold (_type_): _description_
            apply_on (str, optional): _description_. Defaults to 'cal'.
            dimension (list, optional): _description_. Defaults to ['time','channel_idx'].

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        if apply_on=='cal':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.date, spectro)
                self.calibrated_data[spectro] = GROMORA_library.find_bad_channels_stdTb(self.calibrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.calibrated_data
        elif apply_on=='int':
            for spectro in spectrometers:
                bad_channels = self.return_bad_channels(self.date, spectro)
                self.integrated_data[spectro] = GROMORA_library.find_bad_channels_stdTb(self.integrated_data[spectro], bad_channels, stdTb_threshold, dimension)
            return self.integrated_data
        else:
            raise ValueError()
    
    def plot_level2_from_tropospheric_corrected_spectra(self, ac, spectro_dataset, retrieval_param, title, figure_list):
        return GROMORA_library.plot_level2_from_tropospheric_corrected_spectra(spectro_dataset, ac, retrieval_param, title, figure_list)
    
    def plot_level2(self, ac, spectro_dataset, retrieval_param, title, figure_list):
        return GROMORA_library.plot_level2(spectro_dataset, ac, retrieval_param, title, figure_list)

    def plot_level2_test_retrieval(self, ac, retrieval_param, title):
        return GROMORA_library.plot_level2_test_retrieval(ac, retrieval_param, title)