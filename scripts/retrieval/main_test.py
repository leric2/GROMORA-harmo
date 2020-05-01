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

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

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
        """ Get the temperature form a specific temperature channel.
        """
        return self._data[f"Temperature.{channel}"]

    def get_hot_load_temperature(self, time):
        """ Get hot load temperature for a specific time.
        """
        raise NotImplementedError("abstract base class")
    
    
    def read_level1b(self):
        print('reading : ', self._filename)
        return retrieval_module.read_level1b(self._filename)
    
    def plot_level1b_TB(self, LVL1B, calibrationCycle):
        plt.plot(LVL1B.frequencies.values,LVL1B.Tb[calibrationCycle].values)
        pass
    
    def plot_meteo_lvl1b(self,METEO):
        return retrieval_module.plot_meteo_level1b(METEO)
    
    def retrieve(self, LVL1B):
        return retrieval_module.retrieve(LVL1B)
                       
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


def run_retrieval(instrument,retrievalParam):
    """
    In this function we call the retrieval process step-by-step
    
    instrument is actually a concrete implementation of the the abstract 
    DataRetrieval class --> will be initiated here in the future
    
    In the end, we will store this function elsewhere along with the classes definition.
    
    """
    
    LVL1B,METEO,globalAttrs=instrument.read_level1b()
    
    instrument.plot_level1b_TB(LVL1B,retrievalParam['calibrationCycle'])
    
    LVL2 = instrument.retrieve(LVL1B)
    
    # OR
    
    LVL2 = DataRetrievalGROSOM.retrieve(instrument, LVL1B)
    
    if retrievalParam["plot_meteo"]:
        instrument.plot_meteo_lvl1b(METEO)
    
    
    
    
    return LVL1B,METEO,globalAttrs


if __name__=="__main__":
    
    # For testing
    basename="/home/eric/Documents/PhD/GROSOM/Level1/"

    instrumentName="GROMOS"

    retrievalParam=dict()
    retrievalParam["calibrationCycle"]=2
    retrievalParam["plot_meteo"] = True

# Reading function
#retrievalTool["reading_level1b"]=reading_level1b.read_level1b

    if instrumentName=="GROMOS":
       filename = basename+"GROMOS_level1b_AC240_2019_02_12"
       instrument = GROMOS_LvL2(filename)
    else:
       filename = basename+"SOMORA_level1b_AC240_2019_04_16"
       instrument = SOMORA_LvL2(filename)
    
    # Check the structure of the file and maybe use it ?
    print(netCDF4.Dataset(filename+".nc").groups.keys())
    
    LVL1B, METEO, globalAttrsLvl1b = run_retrieval(instrument,retrievalParam)
    
    if globalAttrsLvl1b['instument'] != instrumentName:
        print('wrong instrument')
        
    
    #print(LVL1B)
    
# Check if this is the right instument
"""
if attributes["title"] != retrievalTool["instrument"]:
    raise ValueError("The provided instrument does not correspond to the one "
                    "provided !")
"""
    
    