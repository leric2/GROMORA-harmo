#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:58:02 2020

@author: eric

Module with function for the retrieval of Ozone

"""

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

def read_level1b(filenameLevel1b):
    """Example function with types documented in the docstring.
    Description HERE

    Args:
        param1 (int): The first parameter.
        param2 (str): The second parameter.

    Returns:
        bool: The return value. True for success, False otherwise.

    """
    
    DS = xr.open_dataset(
        filenameLevel1b + ".nc",
        group="spectrometer1",
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        #use_cftime=True,
        )

    globalAttributes=xr.open_dataset(filenameLevel1b + ".nc").attrs
    
    METEO=xr.open_dataset(
        filenameLevel1b+".nc",
        group="meteo",
        decode_times=True,
        decode_coords=True,
        )
    
    return DS,METEO,globalAttributes

def plot_meteo_level1b(METEO):
    """Example function with types documented in the docstring.
    Description HERE

    Args:
        param1 (int): The first parameter.
        param2 (str): The second parameter.

    Returns:
        bool: The return value. True for success, False otherwise.

    """
    fig, axs = plt.subplots(3, sharex=True)
    axs[0].plot(METEO.time.values,METEO.air_temperature.values-273.15) 
    axs[0].set_ylabel('$T_{air}$ [°C]')
    axs[0].set_ylim(-10, 30)
    
    axs[1].plot(METEO.time.values,METEO.relative_humidity.values) 
    axs[1].set_ylabel('rel_hum [-]')
    axs[1].set_ylim(0, 1)
    
    axs[2].plot(METEO.time.values,METEO.precipitation.values) 
    axs[2].set_ylabel('rain [mm]')
    axs[2].set_ylim(0, 5)
    
    axs[2].set_xlabel('time [h]')
    
    fig.suptitle("Meteo Data")
    

def retrieve(level1b):
    for i in level1b.time:
        print('retrieving day : ',i.values)
    return i

if __name__=="__main__":
    main()