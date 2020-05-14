#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:28:46 2020

@author: eric

Collection of functions for dealing with GROSOM data

Including : 
    * level1a, level1b and level2 data
    * a-priori data
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
        
    return DS, METEO, globalAttributes

def apply_correction(level1b_dataset, meteo_ds):
    '''
    Function coming from pywdp, in "troposphere.py" and adapted to GROSOM

    Parameters
    ----------
    
    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    '''
    #Ttrop = get_tropospheric_temperature(ds, mean_temperature_func)
    tau = level1b_dataset.trospheric_opacity.values
    T_trop = meteo_ds.air_temperature.values
    
    Tb_tropospheric_corr = np.ones((len(level1b_dataset.time),len(level1b_dataset.Tb[1])))*np.nan
    for i in range(len(level1b_dataset.time)):
        Tb_tropospheric_corr[i,:] = (level1b_dataset.Tb[i].values - T_trop[i] * (1 - xr.ufuncs.exp(-tau[i]))) / xr.ufuncs.exp(-tau[i])
    
    tbco = xr.DataArray(Tb_tropospheric_corr, dims = ['time', 'channel_idx'])
    
    level1b_dataset = level1b_dataset.assign(Tb_trop_corr = tbco)
    
    return level1b_dataset

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