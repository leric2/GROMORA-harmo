#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Main script for GROSOM retrieval

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

from dotenv import load_dotenv

# For ARTS, we need to specify some paths
#load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

import data_GROSOM
import retrieval_module

#%%

class Level1bDataProcessingGROSOM(ABC):
    '''
    create another class for data processing of level1b
    
    '''
    def __init__(self, filename=None):
        self._filename = filename
    
class DataRetrievalGROSOM(ABC):
    def __init__(self, filename_level1b = None, date=None):
        self._filename_lvl1b = filename_level1b
        self.date = date

        # some attributes for the class (_sss) should be used internally
        self._data = dict()
        self._data["Temperature.0"] = 10
        self._data["Temperature.1"] = 20
        self._data["Temperature.2"] = 30

    # methods for this class
    def get_temerature_reading(self, time, channel):
        ''' Get the temperature form a specific temperature channel.
        '''
        return self._data[f"Temperature.{channel}"]

    def get_hot_load_temperature(self, time):
        ''' Get hot load temperature for a specific time.
        '''
        raise NotImplementedError("abstract base class")
    
    def read_level1b(self):
        ''' 
        Reading level1b dataset and completing the information on the instrument
        
        '''

        print('reading : ', self._filename_lvl1b)

        self.level1b_ds, self.flags, self.meteo_ds, global_attrs_level1b = data_GROSOM.read_level1b(self._filename_lvl1b)

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
        self.raw_fillabview_logfile_warninge_warning = global_attrs_level1b['labview_logfile_warning']

        # some information from the ds
        self.number_of_channels = len(self.level1b_ds.channel_idx.values)

        self.frequencies = self.level1b_ds.frequencies.values
        self.IF = self.level1b_ds.intermediate_freq.values

        self.time = self.level1b_ds.time.values
        self.number_of_time_records = len(self.time)

        return self.level1b_ds, self.flags, self.meteo_ds
    
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
    
    def retrieve_cycle(self, retrieval_param):
        ''' 
        Performing single retrieval for a given calibration cycle (defined in retrieval_param) 
        '''
        return retrieval_module.retrieve_cycle(self, retrieval_param)
    
    def retrieve_cycle_tropospheric_corrected(self, retrieval_param, f_bin = None, tb_bin = None):
        ''' 
        Performing single retrieval for a given calibration cycle uncluding a tropospheric correction
        '''

        if f_bin is not None:
            retrieval_param["binned_ch"] = True
        else:
            retrieval_param["binned_ch"] = False

        return retrieval_module.retrieve_cycle_tropospheric_corrected(self, retrieval_param)

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
    
    def plot_level2_from_tropospheric_corrected_spectra(self, ac, retrieval_param, title):
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
        return data_GROSOM.plot_level2_from_tropospheric_corrected(self.level1b_ds, ac, retrieval_param, title)
    
    def plot_level2(self, ac, retrieval_param, title):
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
        return data_GROSOM.plot_level2(self.level1b_ds, ac, retrieval_param, title)

class GROMOS_LvL2(DataRetrievalGROSOM):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, int_time):
        '''
        Some specific parameters to implement for the GROMOS instances (only constant stuff...)
        '''
        self.observation_frequency = 1.4217504e11
        self.bandwidth = 1e9

        self.spectrometer = "AC240"
        self.instrument_name = "GROMOS"
        
        self.level1_folder = os.path.join(basename_lvl1, self.instrument_name)
        self.level2_folder =  os.path.join(basename_lvl2, self.instrument_name)
        
        self.filename_level1a = os.path.join(
            self.level1_folder,
            self.instrument_name + "_level1a_" +
            self.spectrometer + "_" + date.strftime('%Y_%m_%d')
        )

        if int_time == 1:
            self.filename_level1b = os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1b_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )

            self.filename_level2 = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )
        else:
            self.filename_level1b = os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1b_"+ str(int_time) +"h_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )
    
            self.filename_level2 = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" + str(int_time) +"h_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )

        super().__init__(self.filename_level1b, date)
    
    @staticmethod
    def return_bad_channel_GROMOS(date):
        '''
        to get the bad channels as a function of the date for GROMOS
        
        Parameters
        ----------
        date : datetime object
            DESCRIPTION.
        
        '''
        #if year == 2019,....
        bad_channels = np.arange(16350,16420)

        return bad_channels

    def get_hot_load_temperature(self, time):
        """ On GROMOS, hot load temperature is on channel 0.
        """
        return self.get_temerature_reading(time, 0)

    def find_bad_channels(self, Tb_min, Tb_max, boxcar_size, boxcar_thresh):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        bad_channels = self.return_bad_channel_GROMOS(self.date)
        #self.level1b_ds = super().find_bad_channels(bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
        #

        # ATTENTION, if we return directly the function, the self.level1b_ds is not updated ! TOTHINK !
        self.level1b_ds = data_GROSOM.find_bad_channels(self.level1b_ds, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
        return self.level1b_ds
    
class SOMORA_LvL2(DataRetrievalGROSOM):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, int_time):
        '''
        Some specific parameters to implement for the SOMORA instances (only constant stuff...)
        '''
        self.observation_frequency = 1.4217504e11
        self.bandwidth = 1e9

        self.spectrometer = "AC240"
        self.instrument_name = "SOMORA"

        self.level1_folder = os.path.join(basename_lvl1, self.instrument_name)
        self.level2_folder =  os.path.join(basename_lvl2, self.instrument_name)
        
        self.filename_level1a = os.path.join(
            self.level1_folder,
            self.instrument_name + "_level1a_" + self.spectrometer + "_" + date.strftime('%Y_%m_%d')
        )
        
        if int_time == 1:
            self.filename_level1b = os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1b_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )

            self.filename_level2 = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )
        else:
            self.filename_level1b = os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1b_"+ str(int_time) +"h_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )
    
            self.filename_level2 = os.path.join(
                self.level2_folder,
                self.instrument_name + "_level2_" + str(int_time) +"h_" +
                self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )

        super().__init__(self.filename_level1b, date)

class MOPI5_LvL2(DataRetrievalGROSOM):
    '''
    Implementing the Dataretrieval class for the MOPI5 case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    def __init__(self, date, basename_lvl1, basename_lvl2, int_time):
        '''
        Some specific parameters to implement for the SOMORA instances (only constant stuff...)
        '''
        self.observation_frequency = 110.836e9
        bandwidth = [1.6e9, 1e9, 200e6]

        self.instrument_name = "mopi5"

        self.level1_folder = basename_lvl1
        self.level2_folder =  basename_lvl2

        mopi5_instrument = dict()
        for i, s in enumerate(["U5303", "AC240", "USRP-A"]):
            self.spectrometer = s
            self.bandwidth = bandwidth[i]

            self.filename_level1a = os.path.join(
                self.level1_folder,
                self.instrument_name + "_level1a_" + self.spectrometer + "_" + date.strftime('%Y_%m_%d')
            )

            if int_time == 1:
                self.filename_level1b = os.path.join(
                    self.level1_folder,
                    self.instrument_name + "_level1b_" +
                    self.spectrometer + "_" + date.strftime('%Y_%m_%d')
                )

                self.filename_level2 = os.path.join(
                    self.level2_folder,
                    self.instrument_name + "_level2_" +
                    self.spectrometer + "_" + date.strftime('%Y_%m_%d')
                )
            else:
                self.filename_level1b = os.path.join(
                    self.level1_folder,
                    self.instrument_name + "_level1b_"+ str(int_time) +"h_" +
                    self.spectrometer + "_" + date.strftime('%Y_%m_%d')
                )

                self.filename_level2 = os.path.join(
                    self.level2_folder,
                    self.instrument_name + "_level2_" + str(int_time) +"h_" +
                    self.spectrometer + "_" + date.strftime('%Y_%m_%d')
                )

            mopi5_instrument[s] = super().__init__(self.filename_level1b, date)

    def get_hot_load_temperature(self, time):
        """ On SOMORA, hot load temperature is on channel 2.
        """
        return self.get_temerature_reading(time, 2)    
    
    @staticmethod
    def return_bad_channel_SOMORA(date):
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
    
    def find_bad_channels(self, Tb_min, Tb_max, boxcar_size, boxcar_thresh):
        '''
        Parameters
        ----------
        level1b_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        bad_channels = self.return_bad_channel_SOMORA(self.date)
        #self.level1b_ds = super().find_bad_channels(bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
        #
        return data_GROSOM.find_bad_channels(self.level1b_ds, bad_channels, Tb_min, Tb_max, boxcar_size, boxcar_thresh)
    
def run_retrieval(instrument,retrieval_param):
    '''
    In this function we call the retrieval process step-by-step
    
    instrument is actually a concrete implementation of the the abstract 
    DataRetrieval class --> will be initiated here in the future
    
    In the end, we will store this function elsewhere along with the classes definition.
    
    
    
    level1b_dataset,meteo_ds,global_attrs_level1b = instrument.read_level1b()
    
    instrument.plot_level1b_TB(level1b_dataset,retrieval_param['calibrationCycle'])
    
    level2 = instrument.retrieve(level1b_dataset,retrieval_param)
    
    # OR
    #LVL2 = DataRetrievalGROSOM.retrieve(instrument, level1b_dataset)
    
    if retrieval_param["plot_meteo_ds"]:
        instrument.plot_meteo_ds_level1b_dataset(meteo_ds)
    
    
    return level1b_dataset,meteo_ds,global_attrs_level1b
    '''

def save_single_pdf(filename, figures):
    """
    Save all `figures` to a single PDF. from Jonas
    """
    with PdfPages(filename) as pdf:
        for fig in figures:
            pdf.savefig(fig)


if __name__=="__main__":
    
    instrument_name="mopi5"
    date = datetime.date(2019,2,4)
    int_time = 24

    #basename="/home/eric/Documents/PhD/GROSOM/Level1/"
    #level2_data_folder = "/home/eric/Documents/PhD/GROSOM/Level2/"
    
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/"
    
    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
    #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"

    if instrument_name=="GROMOS":
        instrument = GROMOS_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    elif instrument_name=="SOMORA":
        instrument = SOMORA_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    elif instrument_name=="mopi5":
        basename_lvl1 = "/scratch/MOPI5/Level1/"
        basename_lvl2 = "/scratch/MOPI5/Level2/"
        instrument = MOPI5_LvL2(date, basename_lvl1, basename_lvl2, int_time)
    
    # Dictionnary containing all EXTERNAL retrieval parameters 
    retrieval_param = dict()

    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    retrieval_param["retrieval_type"] = 1

    retrieval_param["obs_freq"] = instrument.observation_frequency
    
    retrieval_param["plot_meteo_ds"] = True
    retrieval_param["number_of_freq_points"] = 601
    retrieval_param["show_f_grid"] = True

    retrieval_param["z_top_sim_grid"] = 97e3
    retrieval_param["z_bottom_sim_grid"] = 800
    retrieval_param["z_resolution_sim_grid"] = 1e3

    retrieval_param["z_top_ret_grid"] = 95e3
    retrieval_param["z_bottom_ret_grid"] = 800
    retrieval_param["z_resolution_ret_grid"] = 3e3

    retrieval_param['increased_var_factor'] = 10
    retrieval_param['unit_var_y']  = 3

    retrieval_param["surface_altitude"] = 10e3

    retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
    retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    retrieval_param["apriori_O3_cov"] = 1.5e-6

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
    
    level1b_ds, flag_ds, meteo_ds = instrument.read_level1b()
    assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'
    
    #retrieval_param = {**global_attrs_level1b, **retrieval_param}
    #else :
    #    raise ValueError('incoherent instrument definition')
    
    #level1b_dataset = instrument.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    
    instrument.level1b_ds = instrument.find_bad_channels(
        Tb_min=0, 
        Tb_max=260, 
        boxcar_size=128, 
        boxcar_thresh=7
    )
    
    # instrument.plot_meteo_ds_level1b_dataset()

    #level1b_dataset = instrument.smooth_corr_spectra(level1b_dataset, retrieval_param)
    #f_sim, y_sim = instrument.forward_model(retrieval_param)
    #plt.plot(f_sim, y_sim[0], level1b_dataset.frequencies.values, level1b_dataset.Tb[1].values)
    
    retrieval_param["integration_cycle"] = 1

    #bin_vector = instrument.create_binning(
    #    instrument.frequencies,
    #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
    #    retrieval_param
    #)

    #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)


    if retrieval_param["retrieval_type"] == 1:
        retrieval_param["observation_altitude"] =  15e3
        ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(retrieval_param, f_bin=None, tb_bin=None)
        figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(ac, retrieval_param, title = 'retrieval_trop_corr')
        level2 = ac.get_level2_xarray()
        level2.to_netcdf(path = instrument.filename_level2+'_'+str(retrieval_param["integration_cycle"])+'.nc')
        save_single_pdf(instrument.filename_level2+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
    elif retrieval_param["retrieval_type"] == 2:
        retrieval_param["observation_altitude"] =  1e3
        ac, retrieval_param = instrument.retrieve_cycle(retrieval_param)
        figure_list = instrument.plot_level2(ac, retrieval_param, title = 'retrieval_h2o')
        save_single_pdf(instrument.filename_level2+'_'+str(retrieval_param["integration_cycle"])+'Perrin_with_h2o.pdf', figure_list)
    else:
        pass
