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
from dotenv import load_dotenv

import GROSOM_library

# %%
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.4')
load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')
# ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
# ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
# ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

# if __name__ == "__main__":

instrument_name = "SOMORA"

# date = pd.date_range(start='2019-01-03', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

# date = pd.date_range(start='2019-01-30', end='2019-06-18')

date = pd.date_range(start='2019-02-21', end='2019-02-21')
#date = datetime.date(2016,1,2)
#date = [datetime.date(2019,3,11), datetime.date(2019,4,3)]

int_time = 1

plot_cycle = False
df_bins = 200e3

plot_all = False
plot_o3_ts = False
plot_selected = False

compare_MERRA2 = True
compare_ECMWF = True

integration_strategy = 'classic'
classic = np.arange(1, 24)

cycle = 14

# %%

basename_lvl1 = "/scratch/GROSOM/Level1/"
basename_lvl2 = "/scratch/GROSOM/Level2/GROMOS/"
if instrument_name=="GROMOS":
    import gromos_classes as gc

    instrument = gc.GROMOS_LvL2(
        date = date, 
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name=="SOMORA":
    import somora_classes as sm
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/SOMORA/"
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
    # o3_test = xr.DataArray(
    #     data=ozone.o3_x.data,
    #     dims={'time','altitude'},
    #     coords=dict(
    #         time=ozone.time.values,
    #         altitude=(['time','altitude'],ozone.o3_z.values)
    #     )
    # )
    o3 = o3.assign_coords({'altitude':ozone.o3_z})
    o3_hourly = o3.resample(time='H', skipna=True).nearest(tolerance='1H')
    o3_hourly.coords['o3_p'] = o3_hourly.coords['o3_p']/100
    o3_hourly.data = o3_hourly.data*1e6
    #o3=o3.swap_dims({'o3_p':'altitude'})
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
    #o3.plot(x='time', y='altitude')
    o3_hourly.plot(
        x='time',
        vmin=0,
        vmax=8,
        cmap='viridis',
        cbar_kwargs={"label": "ozone [PPM]"}
    )
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylabel('P [hPa]')
    plt.tight_layout()
    #o3.plot.imshow(x='time')
    fig.savefig(instrument.level2_folder+'/'+instrument.basename_plot_level2+'ozone_ts_16.pdf')

if plot_selected:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2+instrument.datestr +'_plot_sel_test'
    GROSOM_library.plot_O3_all(
        level2_dataset,
        outname,
        cycles=[0]
    )


if compare_MERRA2:
    merra2_basename = '/mnt/tub/atmosphere/MERRA2/BRN/'
    filename_merra2 =[
        os.path.join(merra2_basename,'MERRA2_BRN_2016_01_diagnostic.h5'),
        os.path.join(merra2_basename,'MERRA2_BRN_2016_02_diagnostic.h5'),
        os.path.join(merra2_basename,'MERRA2_BRN_2016_03_diagnostic.h5')
    ] 
    
    o3_merra2_tot = xr.Dataset()
    counter = 0
    for f in filename_merra2:
        merra2_info = xr.open_dataset(
            f,
            group='info',
            decode_times=False,
            decode_coords=True,
            #use_cftime=True,
        )
        #merra2_info.time.attrs = attrs
        merra2_info['time'] = merra2_info.time - merra2_info.time[0]
    #   # construct time vector

        time_merra2=[]
        for i in range(len(merra2_info.time)):
            time_merra2.append(
                datetime.datetime(
                    int(merra2_info.isel(phony_dim_1=0).year.data[i]), 
                    int(merra2_info.isel(phony_dim_1=0).month.data[i]),  
                    int(merra2_info.isel(phony_dim_1=0).day.data[i]), 
                    int(merra2_info.isel(phony_dim_1=0).hour.data[i]), 
                    int(merra2_info.isel(phony_dim_1=0)['min'].data[i]), 
                    int(merra2_info.isel(phony_dim_1=0).sec.data[i])
                )
            )
        merra2_info['datetime'] = time_merra2
        merra2_decoded = xr.decode_cf(merra2_info)

        merra2 = xr.open_dataset(
            f,
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
        if counter ==0:
            o3_merra2_tot = o3_merra2
        else:
            o3_merra2_tot = xr.concat([o3_merra2_tot, o3_merra2], dim='datetime')
        counter = counter + 1
        
    fig,ax = plt.subplots(1,1)
    o3_merra2_tot.plot(x='datetime', y='altitude', ax=ax,vmin=0, vmax=15)
    ax.set_ylim(5,65)
    #o3_merra2.assign_coords({'altitude':merra2_info.alt.isel(phony_dim_1=0)})
    plt.tight_layout()
    fig.savefig(instrument.level2_folder+'/ozone_ts_16_merra2.pdf')
# %%

if compare_ECMWF:
    ECMWF_folder = '/home/eric/Documents/PhD/ECMWF/'
    counter = 0
    for d in date:
        ECMWF_file = os.path.join(ECMWF_folder, 'ecmwf_oper_v2_BERN_'+d.strftime('%Y%m%d')+'.nc')

        ecmwf = xr.open_dataset(
            ECMWF_file,
            decode_times=True,
            decode_coords=True,
            use_cftime=False,
        )

        if counter == 0:
            ecmwf_ts = ecmwf
        else:
            ecmwf_ts = xr.concat([ecmwf_ts, ecmwf], dim='time')

        counter = counter + 1

    o3_ecmwf = ecmwf_ts.isel(loc=0).ozone_mass_mixing_ratio
    o3_ecmwf['pressure'] = ecmwf_ts['pressure'].isel(loc=0)
    o3_ecmwf.data = o3_ecmwf.data * 1e6
    
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
    #o3.plot(x='time', y='altitude')
    o3_ecmwf.plot(
        x='time',
        y='pressure',
        vmin=0,
        vmax=12,
        cmap='viridis',
        cbar_kwargs={"label": "ozone [PPM]"}
    )
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylabel('P [hPa]')
    plt.tight_layout()
    #o3.plot.imshow(x='time')
    fig.savefig(instrument.level2_folder+'/'+instrument.basename_plot_level2+'ozone_ts_19_ecmwf.pdf')