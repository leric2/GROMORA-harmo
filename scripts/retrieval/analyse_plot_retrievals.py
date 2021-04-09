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
import scipy.io
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from utils_GROSOM import save_single_pdf

from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.lines import Line2D
from dotenv import load_dotenv
#import colormaps as cmaps
import GROSOM_library

#from cmcrameri import cm

# %%
load_dotenv('/home/es19m597/Documents/ARTS/.env.birg-arts24')

# ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
# ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
# ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']
# from apriori_data_GROSOM import read_add_geopotential_altitude
# if __name__ == "__main__":

instrument_name = "SOMORA"

# date = pd.date_range(start='2019-01-03', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

# date = pd.date_range(start='2019-01-30', end='2019-06-18')

date = pd.date_range(start='2017-03-16', end='2017-08-31')
#date = datetime.date(2016,1,2)
#date = [datetime.date(2019,3,11), datetime.date(2019,4,3)]

int_time = 1

plot_cycle = False
df_bins = 200e3

plot_all = False
plot_o3_ts = False
plot_selected = True

compare = False

compare_MERRA2 = False
compare_ECMWF = False
compare_MLS = False

integration_strategy = 'classic'
classic = np.arange(1, 24)

cycle = 14

# %%

colormap = 'viridis'

def read_mls(d1, d2):
    MLS_basename = '/home/esauvageat/Documents/AuraMLS/'
    filename_MLS = os.path.join(MLS_basename,'aura-ozone-at-Bern.mat')        
    
    mls_data = scipy.io.loadmat(filename_MLS)

    ozone_mls = mls_data['o3']
    p_mls = mls_data['p'][0]
    time_mls = mls_data['tm'][0]
    datetime_mls=[]
    for i, t in enumerate(time_mls):
        datetime_mls.append(datetime.datetime.fromordinal(int(t)) + datetime.timedelta(days=time_mls[i]%1) - datetime.timedelta(days = 366))


    ds_mls = xr.Dataset(
        data_vars=dict(
            o3 = (['time','p'], ozone_mls)
        ),
        coords=dict(
            lon=mls_data['longitude2'][0],
            lat= mls_data['latitude2'][0],
            time = datetime_mls,
            p = p_mls
        ),
        attrs=dict(description='ozone time series at bern')
    )

    return ds_mls.sel(time=slice(d1, d2))

if instrument_name=="GROMOS":
    import gromos_classes as gc
    basename_lvl1 = "/scratch/GROSOM/Level1/GROMOS/"    
    basename_lvl2 = "/home/es19m597/Documents/GROMORA/Data/"
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
    basename_lvl2 = "/home/es19m597/Documents/GROMORA/Data/"
    instrument = sm.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name=="compare":
    import somora_classes as sm
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    somora = sm.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
    import gromos_classes as gc
    basename_lvl1 = "/scratch/GROSOM/Level1/GROMOS/"    
    basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    gromos = gc.GROMOS_LvL2(
        date = date, 
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )

if instrument_name=="compare":
    level2_somora = somora.read_level2(
        spectrometers = ['AC240'],
        extra_base='_23'
        )
    level2_gromos = gromos.read_level2(
        spectrometers = ['AC240'],
        extra_base='_23'
        )
    F0 =  somora.observation_frequency
else:
    # Plotting part
    level2_dataset = instrument.read_level2(
        spectrometers = ['AC240'],
        extra_base='_23'
        )
    F0 =  instrument.observation_frequency

if compare:
    ozone_somora = level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_x*1e6
    ozone_gromos = level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_x*1e6
    gromos_error = level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_eo*1e6# +level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_es*1e6
    somora_error = level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_eo*1e6# +level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_es*1e6

    ozone_somora.coords['o3_p'] = ozone_somora.coords['o3_p']/100
    ozone_gromos.coords['o3_p'] = ozone_gromos.coords['o3_p']/100
    gromos_error.coords['o3_p'] = gromos_error.coords['o3_p']/100
    somora_error.coords['o3_p'] = somora_error.coords['o3_p']/100

    mls = read_mls(date[0].strftime('%Y-%m-%d'), date[-1].strftime('%Y-%m-%d'))
    o3_mls = mls.o3
    rel_diff = 100*(ozone_somora.mean(dim='time') - ozone_gromos.mean(dim='time'))/ozone_somora.mean(dim='time')
   # rel_diff_gromos_mls = 100*(o3_mls.mean(dim='time') - ozone_gromos.mean(dim='time'))/o3_mls.mean(dim='time')
   # rel_diff_somora_mls = 100*(o3_mls.mean(dim='time') - ozone_somora.mean(dim='time'))/o3_mls.mean(dim='time')

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(9,6))

    ozone_somora.mean(dim='time').plot(y='o3_p',ax=axs[0],color='red', label='SOMORA')
    ozone_gromos.mean(dim='time').plot(y='o3_p',ax=axs[0],color='blue', label='GROMOS')
    o3_mls.mean(dim='time').plot(y='p',ax=axs[0],color='black', label='MLS')

    gromos_plus = ozone_gromos.mean(dim='time') + gromos_error.mean(dim='time')
    gromos_minus = ozone_gromos.mean(dim='time') - gromos_error.mean(dim='time')
    somora_plus = ozone_somora.mean(dim='time') + somora_error.mean(dim='time')
    somora_minus = ozone_somora.mean(dim='time') - somora_error.mean(dim='time')

    gromos_minus.plot(y='o3_p',ax=axs[0], ls='--', color='blue', alpha=0.4, label='mean obs error GROMOS')
    gromos_plus.plot(y='o3_p',ax=axs[0], ls='--', color='blue', alpha=0.4)

    somora_plus.plot(y='o3_p',ax=axs[0], ls='--', color='red', alpha=0.4, label='mean obs error, SOMORA')
    somora_minus.plot(y='o3_p',ax=axs[0], ls='--', color='red', alpha=0.4)

   # apriori_mean.mean(dim='time').plot(y='o3_p',label='GROMOS')
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].legend()
    axs[0].set_ylim(500,1e-3)
    axs[0].set_ylabel('pressure [hPa]')
    axs[0].grid(axis='x', linewidth=0.5)

    rel_diff.plot(y='o3_p',ax=axs[1], color='green', label='(SOM-GRO)/SOM')
  #  rel_diff_gromos_mls.plot(y='o3_p',ax=axs[1], color='blue', label='(MLS-GRO)/MLS')
  #  rel_diff_somora_mls.plot(y='o3_p',ax=axs[1], color='red', label='(MLS-SOM)/MLS')
    axs[1].set_xlim(-75,75)
    axs[1].set_xlabel('relative difference [%]')
    axs[1].legend() 
    axs[1].grid(axis='x', linewidth=0.5)

    for a in axs:
        a.grid(which='both', axis='y', linewidth=0.5)
    plt.suptitle('Ozone comparison with '+ str(len(date)) +' days ' + pd.to_datetime(ozone_somora.time.mean().data).strftime('%Y-%m-%d %H:%M'))
    plt.tight_layout()

    # save mean profiles for these dates
    o3_mls.mean(dim='time').to_netcdf(somora.level2_folder+'/'+'mls_mean_o3_'+date.mean().strftime('%Y-%m-%d')+'.nc')
    ozone_gromos.mean(dim='time').to_netcdf(somora.level2_folder+'/'+'gromos_mean_o3_'+date.mean().strftime('%Y-%m-%d')+'.nc')
    ozone_somora.mean(dim='time').to_netcdf(somora.level2_folder+'/'+'somora_mean_o3_'+date.mean().strftime('%Y-%m-%d')+'.nc')

    fig.savefig(somora.level2_folder+'/'+'ozone_comparison_'+pd.to_datetime(ozone_somora.time.mean().data).strftime('%Y-%m-%d')+'.pdf')

if plot_all:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2+instrument.datestr +'_plot_all_test_polyfit2'
    GROSOM_library.plot_O3_all(level2_dataset, outname)

if plot_o3_ts:
    ozone = level2_dataset['AC240'].isel(o3_lat=0, o3_lon=0)
    o3 = ozone.o3_x
    o3['altitude'] =  ozone.o3_z / 1e3
    # o3_test = xr.DataArray(
    #     data=ozone.o3_x.data,
    #     dims={'time','altitude'},
    #     coords=dict(
    #         time=ozone.time.values,
    #         altitude=(['time','altitude'],ozone.o3_z.values)
    #     )
    # )
   # o3 = o3.assign_coords({'altitude':ozone.o3_z})
    o3_hourly = o3.resample(time='H', skipna=True).nearest(tolerance='1H')
    o3.coords['o3_p'] = o3.coords['o3_p']/100
    o3.data = o3.data*1e6
   # o3=o3.swap_dims({'o3_p':'altitude'})
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
  #  o3_hourly = o3_hourly.swap_dims({'o3_p':'geometric_height'})
   # o3_hourly['geometric_height'] = o3_hourly.o3_z
    #o3.plot(x='time', y='altitude')
    o3.plot(
        x='time',
        y='o3_p',
        vmin=0,
        vmax=9,
        cmap=colormap,
        cbar_kwargs={"label": "ozone [PPM]"}
    )
    ax.set_ylim(0.01,500)
    ax.invert_yaxis()
    ax.set_yscale('log')
    
    ax.set_ylabel('Pressure [hPa]')
    plt.tight_layout()
    #o3.plot.imshow(x='time')
    fig.savefig(instrument.level2_folder+'/'+instrument.basename_plot_level2+'ozone_ts_01-05-2019.pdf')

if compare_MLS:
    MLS_basename = '/home/esauvageat/Documents/AuraMLS/'
    filename_MLS = os.path.join(MLS_basename,'aura-ozone-at-Bern.mat')        
    
    mls_data = scipy.io.loadmat(filename_MLS)

    ozone_mls = mls_data['o3']
    p_mls = mls_data['p'][0]
    time_mls = mls_data['tm'][0]
    datetime_mls=[]
    for i, t in enumerate(time_mls):
        datetime_mls.append(datetime.datetime.fromordinal(int(t)) + datetime.timedelta(days=time_mls[i]%1) - datetime.timedelta(days = 366))


    ds_mls = xr.Dataset(
        data_vars=dict(
            o3 = (['time','p'], ozone_mls)
        ),
        coords=dict(
            lon=mls_data['longitude2'][0],
            lat= mls_data['latitude2'][0],
            time = datetime_mls,
            p = p_mls
        ),
        attrs=dict(description='ozone time series at bern')
    )
    #ds_mls = xr.decode_cf(ds_mls)
    #ds_mls.time.encoding['units'] = 'seconds since 1970-01-01 00:00:00'
    #ds_mls.time.encoding['calendar'] = "proleptic_gregorian"

    #ds_mls.to_netcdf('/home/esauvageat/Documents/AuraMLS/ozone_bern_ts.nc', format='NETCDF4')

    o3_mls = ds_mls.sel(time=slice("2019-01-10", "2019-05-26")).o3
    
    monthly_mls = o3_mls.resample(time='1m', skipna=True).mean()

    fig,ax = plt.subplots(1,1)
    #monthly_mls.plot(x='time', y='p', ax=ax ,vmin=0, vmax=9)
    o3_mls.plot(x='time', y='p', ax=ax ,vmin=0, vmax=9, cmap=colormap)
    ax.set_yscale('log')
    ax.set_ylim(0.01,500)
    ax.invert_yaxis()
    ax.set_ylabel('P [hPa]')
    plt.tight_layout()
    fig.savefig(instrument.level2_folder+'/ozone_mls_01-04-2019.pdf')


if plot_selected:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2+instrument.datestr +'_plot_sel_polyfit2'
    GROSOM_library.plot_O3_all(
        level2_dataset,
        outname,
        cycles=[1,3,6,9]
    )


if compare_MERRA2:
    merra2_basename = '/storage/tub/atmosphere/MERRA2/BRN/'
    filename_merra2 =[
       # os.path.join(merra2_basename,'MERRA2_BRN_2016_01_diagnostic.h5'),
        os.path.join(merra2_basename,'MERRA2_BRN_2016_02_diagnostic.h5'),
        os.path.join(merra2_basename,'MERRA2_BRN_2016_03_diagnostic.h5'),
        os.path.join(merra2_basename,'MERRA2_BRN_2016_04_diagnostic.h5')
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
            group='trace_gas',        
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

    o3_merra2_tot = o3_merra2_tot.sel(datetime=slice("2016-02-15", "2016-04-07"))
        
    fig,ax = plt.subplots(1,1)
    o3_merra2_tot.plot(x='datetime', y='altitude', ax=ax,vmin=0, vmax=15)
    ax.set_ylim(5,65)
    #o3_merra2.assign_coords({'altitude':merra2_info.alt.isel(phony_dim_1=0)})
    plt.tight_layout()
    fig.savefig(instrument.level2_folder+'/ozone_ts_16_merra2.pdf')
# %%

if compare_ECMWF:
    ECMWF_folder = '/scratch/ECMWF/'
    counter = 0
    for d in date:
        ECMWF_file = os.path.join(ECMWF_folder, 'ecmwf_oper_v2_PAYERNE_'+d.strftime('%Y%m%d')+'.nc')

        ecmwf_og = xr.open_dataset(
            ECMWF_file,
            decode_times=True,
            decode_coords=True,
            use_cftime=False,
        )
       # ecmwf_og.swap_dims({'level':'pressure'} )
        # for i in range(len(ecmwf_og.time.data)):
        #     ecmwf = ecmwf_og.isel(loc=0, time=i)
        #     ecmwf = read_add_geopotential_altitude(ecmwf)
        if counter == 0:
            ecmwf_ts = ecmwf_og
        else:
            ecmwf_ts = xr.concat([ecmwf_ts, ecmwf_og], dim='time')
        
        counter = counter + 1

    o3_ecmwf = ecmwf_ts.isel(loc=0).ozone_mass_mixing_ratio
    o3_ecmwf['pressure'] = ecmwf_ts['pressure'].isel(loc=0)/100
    o3_ecmwf.data = o3_ecmwf.data * 1e6
    
    fig2 = plt.figure(num=1)
    ax = fig2.subplots(1)
    #o3.plot(x='time', y='altitude')
    o3_ecmwf.plot(
        x='time',
        y='pressure',
        vmin=0,
        vmax=15,
        cmap='viridis',
        cbar_kwargs={"label": "ozone [PPM]"}
    )
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylabel('P [hPa]')
    plt.tight_layout()
    #o3.plot.imshow(x='time')
    fig2.savefig(instrument.level2_folder+'/'+'ozone_ts_16_ecmwf_payerne.pdf')
# %%
    from retrievals.data.ecmwf import ECMWFLocationFileStore
    from retrievals.data.ecmwf import levels
    
    ecmwf_prefix = 'ecmwf_oper_v2_BERN_%Y%m%d.nc'
    t1 = date[0]
    t2 = date[2]
    ecmwf_store = ECMWFLocationFileStore(ECMWF_folder, ecmwf_prefix)
    ds_ecmwf = (
        ecmwf_store.select_time(t1, t2, combine='by_coords')
        .mean(dim='time')
        .swap_dims({"level":"pressure"})
    )

    ds_ecmwf = read_add_geopotential_altitude(ds_ecmwf)
    