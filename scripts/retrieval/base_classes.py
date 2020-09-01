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

import data_GROSOM
import retrieval_module
import mopi5_retrievals

#%%

class MicrowaveRadiometerIAP(ABC):
    '''
    The base abstract class for every IAP MW radiometer. 

    It has multiple children classes, for integration, retrieval, etc...

    '''
    def __init__(
        self,
        instrument_name=None,
        observation_frequency=None,
        spectrometers=None,
        ):
        
        self.instrument_name = instrument_name
        self.observation_frequency = observation_frequency
        self.spectrometers = spectrometers

class Integration(MicrowaveRadiometerIAP):
    '''
    create another class for integration of level1a data
    
    '''
    def __init__(
        self,
        instrument_name=None,
        observation_frequency=None,
        spectrometers=None,
        integration_strategy='simple',
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

        super().__init__(instrument_name, observation_frequency, spectrometers)

        self.dates = dates

        if len(self.dates) > 1:
            self.multiday = True
        else:
            self.multiday = False

        self.calibrated_data = dict()
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

            if self.integration_strategy == 'simple':
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
                    self.instrument_name + "_level1b_"+ self.integration_strategy + '_' + str(self.int_time) +"h_" +
                    s + "_" + self.datestr
                    )

    def read_level1a(self):
        ''' 
        Reading level1a in netCDF format
        
        '''
        if self.multiday:
            raise NotImplementedError('TODO: Implement mutli days level1a reading')
            # Append xarray dataset together ?
        else:
            for s in self.spectrometers:
                print('reading : ', self.filename_level1a[s][0])
                #self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = data_GROSOM.read_level1(self._filename_lvl1b)
                self.calibrated_data[s], self.calibration_flags[s], meteo_data, global_attrs_level1a = data_GROSOM.read_level1(self.filename_level1a[s][0])
        
            self.meteo_complete = meteo_data

        # Meta data
        self.institution = global_attrs_level1a['institution']
        self.instrument_name_from_level1a = global_attrs_level1a['instrument']
        self.location = global_attrs_level1a['location']

        self.number_of_spectrometer = global_attrs_level1a['number_of_spectrometer']
        self.calibration_version = global_attrs_level1a['calibration_version']

        self.raw_data_filename = global_attrs_level1a['raw_data']
        self.raw_data_software_version = global_attrs_level1a['raw_data_software_version']
        
        self.filename_level1a = global_attrs_level1a['filename']
        self.raw_file_warning = global_attrs_level1a['raw_file_warning']
        self.labview_logfile_warning = global_attrs_level1a['labview_logfile_warning']

        return self.calibrated_data, self.calibration_flags, self.meteo_complete

    def integrate(self, spectrometers):
        '''
        Generic function for the integration of Integration class
        '''
        if self.integration_strategy == 'simple':
            print('Performing simple integration')
            raise NotImplementedError()
        elif self.integration_strategy == 'TOD':
            print('Performing integration based on the time of day')
            raise NotImplementedError()
        elif self.integration_strategy == 'MeanTb':
            print('Performing integration based on mean Tb')
            raise NotImplementedError()        
        elif self.integration_strategy == '...':
            print('Performing integration based on ...')
            raise NotImplementedError()

class DataRetrieval(MicrowaveRadiometerIAP):
    def __init__(
        self,
        instrument_name=None,
        observation_frequency=None,
        spectrometers=None,
        integration_time=None,
        date=None,
        level1_folder=None,
        level2_folder=None
        ):
        
        self.instrument_name = instrument_name
        self.date = date
        self.datestr = self.date.strftime('%Y_%m_%d')
        self.observation_frequency = observation_frequency
        self.spectrometers = spectrometers
        self.int_time = integration_time
        self.level1_folder = level1_folder
        self.level2_folder = level2_folder

        super().__init__(instrument_name, observation_frequency, spectrometers)
        # some attributes for the class (_sss) should be used internally
        #self._data = dict()
        #self._data["Temperature.0"] = 10
        #self._data["Temperature.1"] = 20
        #self._data["Temperature.2"] = 30

    # methods for this class
    #def get_temerature_reading(self, time, channel):
    #    ''' Get the temperature form a specific temperature channel.
    #   '''
    #    return self._data[f"Temperature.{channel}"]

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
        self.filename_level1a = dict()
        self.filename_level1b = dict()
        self.filename_level2 = dict()
        for i, s in enumerate(self.spectrometers):
            self.filename_level1a[s] = os.path.join(
            self.level1_folder,
            self.instrument_name + "_level1a_" +
            s + "_" + self.datestr
            )
            if self.int_time == 1:
                self.filename_level1b[s] = os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1b_" +
                s + "_" + self.datestr
                )

                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                s + "_" + self.datestr
                )
            else:
                self.filename_level1b[s] = os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1b_"+ str(self.int_time) +"h_" +
                s + "_" + self.datestr
                )

                self.filename_level2[s] = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" + str(self.int_time) +"h_" +
                s + "_" + self.datestr
                )
        
            print('reading : ', self.filename_level1b[s])
            #self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = data_GROSOM.read_level1(self._filename_lvl1b)
            self.data[s], self.flags[s], meteo_data, global_attrs_level1b = data_GROSOM.read_level1(self.filename_level1b[s])
        
        self.meteo = meteo_data

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
        data_GROSOM.plot_meteo_level1b(self.meteo_ds)
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
        return data_GROSOM.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    
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

        return data_GROSOM.create_bin_vector(self.observation_frequency, freq, tb, n_f, bw_extra)

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

        return data_GROSOM.bin_spectrum(freq, tb, bin_vect)

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
        return data_GROSOM.smooth_corr_spectra(level1b_dataset, retrieval_param)
    
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
        return data_GROSOM.find_bad_channels(self.level1b_ds, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)

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
        self.data[spectro] = data_GROSOM.find_bad_channels_stdTb(self.data[spectro], bad_channels, stdTb_threshold)

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
        return data_GROSOM.plot_level2_from_tropospheric_corrected(spectro_dataset, ac, retrieval_param, title)
    
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
        return data_GROSOM.plot_level2(spectro_dataset, ac, retrieval_param, title)

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
        return data_GROSOM.plot_level2_test_retrieval(ac, retrieval_param, title)

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
