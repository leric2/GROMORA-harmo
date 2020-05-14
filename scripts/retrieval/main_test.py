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

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from dotenv import load_dotenv

# For ARTS, we need to specify some paths
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']

# Following variables could already be set in ~/.bashrc

ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

import data_GROSOM
import retrieval_module

class DataRetrievalGROSOM(ABC):
    def __init__(self, filename=None):
        self._filename = filename

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
        print('reading : ', self._filename)
    
        # self.DS = xr.open_dataset(
        #     self._filename + ".nc",
        #     group="spectrometer1",
        #     mask_and_scale=True,
        #     decode_times=True,
        #     decode_coords=True,
        #     #use_cftime=True,
        #     )
    
        # self.global_attrs = xr.open_dataset(self._filename + ".nc").attrs
        
        # self.meteo=xr.open_dataset(
        #     self._filename+".nc",
        #     group="meteo",
        #     decode_times=True,
        #     decode_coords=True,
        #     )       
        # return self.DS, self.global_attrs, self.meteo
        return data_GROSOM.read_level1b(self._filename)
    
    def plot_level1b_TB(self, level1b_dataset, calibration_cycle):
        plt.plot(level1b_dataset.frequencies.values,level1b_dataset.Tb[calibration_cycle].values)
        pass
    
    def plot_meteo_ds_level1b_dataset(self, meteo_ds):
        return data_GROSOM.plot_meteo_level1b(meteo_ds)
    
    def forward_model(self, retrieval_param):
        ''' 
        Performing a FM according to the retrieval_param dictionnary as input
        parameters
        
        '''
        return retrieval_module.forward_model(retrieval_param)
    
    def retrieve_cycle(self, level1b_dataset, meteo_ds, retrieval_param):
        ''' 
        Performing single retrieval for a given calibration cycle (defined in retrieval_param) 
        '''
        
        return retrieval_module.retrieve_cycle(level1b_dataset, meteo_ds, retrieval_param)
    
    '''
    def freq_grid(self, retrieval_param):
        n_f = retrieval_param["number_of_freq_points"]# Number of points
        bw = retrieval_param["bandwidth"]  # Bandwidth
        x = np.linspace(-1, 1, n_f)
        f_grid = x ** 3 + x / 10
        f_grid = f_grid * bw / (max(f_grid) - min(f_grid))
        return retrieval_param
    '''     

    def apply_correction(self, level1b_dataset, meteo_ds):   
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
        return data_GROSOM.apply_correction(level1b_dataset, meteo_ds)
    
    def plot_level2(self, level1b_dataset, ac, retrieval_param, title):
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
        return retrieval_module.plot(level1b_dataset, ac, retrieval_param, title)

class GROMOS_LvL2(DataRetrievalGROSOM):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    
    def get_hot_load_temperature(self, time):
        """ On GROMOS, hot load temperature is on channel 0.
        """
        return self.get_temerature_reading(time, 0)
    
    def plot_meteo_ds_level1b_dataset(self, meteo_ds):
        print('overcoming this function')
    
class SOMORA_LvL2(DataRetrievalGROSOM):
    '''
    Implementing the Dataretrieval class for the GROMOS case.
    Hereafter we define some parameters and methods specific to this 
    instrument.
    '''
    
    def get_hot_load_temperature(self, time):
        """ On SOMORA, hot load temperature is on channel 2.
        """
        return self.get_temerature_reading(time, 2)    


def run_retrieval(instrument,retrieval_param):
    '''
    In this function we call the retrieval process step-by-step
    
    instrument is actually a concrete implementation of the the abstract 
    DataRetrieval class --> will be initiated here in the future
    
    In the end, we will store this function elsewhere along with the classes definition.
    
    '''
    
    level1b_dataset,meteo_ds,global_attrs_level1b = instrument.read_level1b()
    
    instrument.plot_level1b_TB(level1b_dataset,retrieval_param['calibrationCycle'])
    
    level2 = instrument.retrieve(level1b_dataset,retrieval_param)
    
    # OR
    #LVL2 = DataRetrievalGROSOM.retrieve(instrument, level1b_dataset)
    
    if retrieval_param["plot_meteo_ds"]:
        instrument.plot_meteo_ds_level1b_dataset(meteo_ds)
    
    
    return level1b_dataset,meteo_ds,global_attrs_level1b


def save_single_pdf(filename, figures):
    """
    Save all `figures` to a single PDF. from Jonas
    """
    with PdfPages(filename) as pdf:
        for fig in figures:
            pdf.savefig(fig)



if __name__=="__main__":
    
    # For testing
    basename="/home/eric/Documents/PhD/GROSOM/Level1/"
    
    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
    
    instrument_name="GROMOS"
    
    # Reading function
    #retrievalTool["reading_level1b"]=reading_level1b.read_level1b

    if instrument_name=="GROMOS":
       filename = basename+"GROMOS_level1b_AC240_2019_02_12"
       instrument = GROMOS_LvL2(filename)
    else:
       filename = basename+"SOMORA_level1b_AC240_2019_04_16"
       instrument = SOMORA_LvL2(filename)
    
    retrieval_param = dict()
    retrieval_param["integration_cycle"] = 2
    retrieval_param["plot_meteo_ds"] = True
    retrieval_param["number_of_freq_points"] = 601
    retrieval_param["altitude"] = 461
    retrieval_param["zenith_angle"]=80
    retrieval_param["azimuth_angle"]=32
    retrieval_param["observation_altitude"] = 10e3
    retrieval_param['obs_freq'] = 1.4217504e11
    retrieval_param['line_file'] = line_file
    
    
    fascod_atmosphere = 'midlatitude-summer'
    retrieval_param['prefix_atm'] = ARTS_DATA_PATH + "/planets/Earth/Fascod/{}/{}".format(fascod_atmosphere,fascod_atmosphere)
    
    # Check the structure of the file and maybe use it ?
    #print(netCDF4.Dataset(filename+".nc").groups.keys())
    
    #level1b_dataset, meteo_ds, global_attr_lvl1b_ds = run_retrieval(instrument,retrievalParam)
    
    level1b_dataset, meteo_ds, global_attrs_level1b = instrument.read_level1b()
    
    if global_attrs_level1b['instument'] == instrument_name:
        # merge attrs level1b in retrievalParam
        retrieval_param = {**global_attrs_level1b, **retrieval_param}
    else :
        raise ValueError('incoherent instrument definition')
    
    level1b_dataset = instrument.apply_correction(level1b_dataset, meteo_ds)
    
    #f_sim, y_sim = instrument.forward_model(retrieval_param)
    #plt.plot(f_sim, y_sim[0], level1b_dataset.frequencies.values, level1b_dataset.Tb[1].values)
    
    ac, retrieval_param = instrument.retrieve_cycle(level1b_dataset, meteo_ds, retrieval_param)
    
    figure_list = instrument.plot_level2(level1b_dataset, ac, retrieval_param, 'first try')
    
    save_single_pdf('firstTryO3retrieval.pdf', figure_list)
# Check if this is the right instument
"""
if attributes["title"] != retrievalTool["instrument"]:
    raise ValueError("The provided instrument does not correspond to the one "
                    "provided !")
"""
    
    