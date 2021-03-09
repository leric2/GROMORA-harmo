#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Integration script for IAP instruments

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

from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.lines import Line2D
import plotly_express as px
from dotenv import load_dotenv

import GROSOM_library

# %%
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.4')
# ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
# ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
# ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

# if __name__ == "__main__":

instrument_name = "SOMORA"

# date = pd.date_range(start='2019-01-03', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

# date = pd.date_range(start='2019-01-30', end='2019-06-18')

date = pd.date_range(start='2016-02-16', end='2016-02-16')
#date = datetime.date(2019,1,15)
#date = [datetime.date(2019,3,11), datetime.date(2019,4,3)]

int_time = 1

plot_cycle = False
df_bins = 200e3

plot_all = False
plot_o3_ts = False
plot_selected = False

compare_MERRA2 = True

integration_strategy = 'classic'
classic = np.arange(1, 24)

cycle = 14
# %%
if instrument_name=="GROMOS":
    import gromos_classes as gc
    basename_lvl1 = "/home/eric/Documents/PhD/GROSOM/Data/"
    basename_lvl2 = "/home/eric/Documents/PhD/GROSOM/Data/"  
    instrument = gc.GROMOS_LvL2(
        date = date, 
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name=="SOMORA":
    import somora_classes as sm
    basename_lvl1 = "/home/eric/Documents/PhD/GROSOM/Data/"
    basename_lvl2 = "/home/eric/Documents/PhD/GROSOM/Data/"  
    instrument = sm.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )


 # Plotting part
level2_dataset = instrument.read_level2(
    spectrometers = ['AC240'],
    extra_base='_23'
    )
F0 =  instrument.observation_frequency

if plot_all:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2+instrument.datestr +'_plot_all_test'
    GROSOM_library.plot_O3_all(level2_dataset, outname)

if plot_o3_ts:
    ozone = level2_dataset['AC240'].isel(o3_lat=0, o3_lon=0)
    o3 = ozone.o3_x
    o3 = o3.assign_coords({'altitude':ozone.o3_z})
    o3_hourly = o3.resample(time='H', skipna=False).nearest(tolerance='1H')
    
    
    #o3=o3.swap_dims({'o3_p':'altitude'})
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
    #o3.plot(x='time', y='altitude')
    o3_hourly.plot(x='time')
    ax.invert_yaxis()
    ax.set_yscale('log')
    #o3.plot.imshow(x='time')


if plot_selected:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2+instrument.datestr +'_plot_sel_test'
    GROSOM_library.plot_O3_all(
        level2_dataset,
        outname,
        cycles=[1,6,12,18]
    )


if compare_MERRA2:
    merra2_basename = '/mnt/tub/atmosphere/MERRA2/BRN/'
    filename_merra2 = os.path.join(merra2_basename,'MERRA2_BRN_'+date[0].strftime('%Y_%m')+'_diagnostic.h5')
    
    attrs = {"units": 'days since '+date[0].strftime('%Y-%m-%d')}
    merra2_info = xr.open_dataset(
        filename_merra2,
        group='info',
        decode_times=False,
        decode_coords=True,
        #use_cftime=True,
    )
    #merra2_info.time.attrs = attrs
    merra2_info['time'] = merra2_info.time - merra2_info.time[0]
    #
    # construct time vector
    time_merra2=[]
    for i in range(len(merra2_info.time)):
        time_merra2.append(
            datetime.datetime(
                merra2_info.year.data[i], 
                merra2_info.month.data[i],  
                merra2_info.day.data[i], 
                merra2_info.hour.data[i], 
                merra2_info['min'].data[i], 
                merra2_info.sec.data[i]
            )
        )
    merra2_info['datetime'] = time_merra2
    merra2_decoded = xr.decode_cf(merra2_info)

    merra2 = xr.open_dataset(
        filename_merra2,
        group='wind',        
        decode_times=True,
        decode_coords=True,
        #use_cftime=True,
    )

    o3_merra2 = merra2.O3

    o3_merra2 = o3_merra2.swap_dims({'phony_dim_6':'altitude'})
    o3_merra2['altitude']=merra2_decoded.alt.isel(phony_dim_1=0).data

    o3_merra2 = o3_merra2.swap_dims({'phony_dim_5':'datetime'})
    o3_merra2['datetime']=merra2_decoded.datetime.data

    o3_merra2.data = o3_merra2.data*1e6

    fig,ax = plt.subplots(1,1)
    o3_merra2.plot(x='datetime', y='altitude', ax=ax,vmin=0, vmax=15)
    ax.set_ylim(5,65)
    #o3_merra2.assign_coords({'altitude':merra2_info.alt.isel(phony_dim_1=0)})

# %%
