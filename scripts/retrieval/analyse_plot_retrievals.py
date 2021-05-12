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
import datetime
import os
from abc import ABC

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import scipy.io
import xarray as xr
from dotenv import load_dotenv
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import (AutoMinorLocator, FormatStrFormatter,
                               MultipleLocator)

import GROSOM_library
import mopi5_library
from utils_GROSOM import save_single_pdf

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]})

load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')
# ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
# ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
# ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']
# from apriori_data_GROSOM import read_add_geopotential_altitude
# if __name__ == "__main__":

instrument_name = "compare"

# date = pd.date_range(start='2019-01-03', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

# date = pd.date_range(start='2019-01-30', end='2019-06-18')

date = pd.date_range(start='2018-07-15', end='2018-08-15')
#date = pd.date_range(start='2017-09-01', end='2018-01-05')
#date = datetime.date(2016,1,2)
#date = [datetime.date(2019,3,11), datetime.date(2019,4,3)]

int_time = 1

plot_cycle = False
df_bins = 200e3

plot_all = False
plot_all_mopi5 = False
plot_o3_ts = False
save_o3 = False
plot_selected = False
plot_fshift = False
plot_cost = False
plot_o3_diff_waccm = False
read_waccm_clim = False
compare = True

compare_MERRA2 = False
compare_ECMWF = False
compare_MLS = False

integration_strategy = 'classic'
classic = np.arange(1, 24)

cycle = 14
spectros = ['U5303','AC240','USRP-A'] #
spectros = ['USRP-A','U5303'] 
spectros = ['AC240'] 

ex = 'fascodunbiased_all'
ex = '_fascod_fix_noise_3'
ex = '_waccm_cov_yearly'
#ex = '_waccm'
# %%

colormap = 'cividis'  # 'viridis' #, batlow_map cmap_crameri cividis


def read_mls(d1, d2):
    MLS_basename = '/home/esauvageat/Documents/AuraMLS/'
    filename_MLS = os.path.join(MLS_basename, 'aura-ozone-at-Bern.mat')

    mls_data = scipy.io.loadmat(filename_MLS)

    ozone_mls = mls_data['o3']
    p_mls = mls_data['p'][0]
    time_mls = mls_data['tm'][0]
    datetime_mls = []
    for i, t in enumerate(time_mls):
        datetime_mls.append(datetime.datetime.fromordinal(
            int(t)) + datetime.timedelta(days=time_mls[i] % 1) - datetime.timedelta(days=366))

    ds_mls = xr.Dataset(
        data_vars=dict(
            o3=(['time', 'p'], ozone_mls)
        ),
        coords=dict(
            lon=mls_data['longitude2'][0],
            lat=mls_data['latitude2'][0],
            time=datetime_mls,
            p=p_mls
        ),
        attrs=dict(description='ozone time series at bern')
    )

    return ds_mls.sel(time=slice(d1, d2))


if instrument_name == "GROMOS":
    import gromos_classes as gc
    basename_lvl1 = "/scratch/GROSOM/Level1/GROMOS/"
    basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    instrument = gc.GROMOS_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name == "SOMORA":
    import somora_classes as sm
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    instrument = sm.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name == "mopi5":
    import mopi5_classes as mp
    basename_lvl1 = "/storage/tub/instruments/mopi5/level1"
    basename_lvl2 = "/storage/tub/instruments/mopi5/level2"
    instrument = mp.MOPI5_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name == "compare":
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
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )

if instrument_name == "compare":
    level2_somora = somora.read_level2(
        spectrometers=['AC240'],
        extra_base='_waccm_cov_yearly'
    )
    level2_gromos = gromos.read_level2(
        spectrometers=['AC240'],
        extra_base='_waccm_cov_yearly'
    )
    F0 = somora.observation_frequency
else:
    # Plotting part
    level2_dataset = instrument.read_level2(
        spectrometers=spectros,
        extra_base=ex
    )
    F0 = instrument.observation_frequency


if plot_cost:
    if instrument_name == 'mopi5':
        for s in spectros:
            fig, axs = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(9, 6))
            end_cost = level2_dataset[s].oem_diagnostics[:, 3]
            axs.plot(end_cost)
            fig.savefig(instrument.level2_folder+'/' +
                        instrument.datestr+s+ex+'_end_cost.pdf', dpi=500)
    else:
        for s in spectros:
            fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(9, 6))
            end_cost = level2_dataset[s].oem_diagnostics[:, 3]
            noise_input = level2_dataset[s].median_noise
            end_cost.plot(marker='.', ax=axs[0])
            noise_input.plot(marker='.',ax=axs[1])
            #axs[1].set_ylim(0,1)
          #  end_cost.plot(ax=axs, ylim=(0.75,8))
            fig.savefig(instrument.level2_folder+'/'+instrument.basename_plot_level2 +
                        instrument.datestr+'_end_cost.pdf', dpi=500)

if compare:

    mr_somora = level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_mr.data
    mr_gromos = level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_mr.data
    mr_mean_gromos = level2_somora['AC240'].isel(
        o3_lat=0, o3_lon=0).o3_mr.mean(dim='time').data

    ozone_somora = level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_x*1e6
    ozone_gromos = level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_x*1e6
    # +level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_es*1e6
    gromos_error = level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_eo*1e6
    # +level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_es*1e6
    somora_error = level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_eo*1e6

    ozone_somora.coords['o3_p'] = ozone_somora.coords['o3_p']/100
    ozone_gromos.coords['o3_p'] = ozone_gromos.coords['o3_p']/100
    gromos_error.coords['o3_p'] = gromos_error.coords['o3_p']/100
    somora_error.coords['o3_p'] = somora_error.coords['o3_p']/100

    mls = read_mls(date[0].strftime('%Y-%m-%d'), date[-1].strftime('%Y-%m-%d'))
    o3_mls = mls.o3

    rel_diff = 100*(ozone_somora.mean(dim='time') -
                    ozone_gromos.mean(dim='time'))/ozone_somora.mean(dim='time')
   # rel_diff_gromos_mls = 100*(o3_mls.mean(dim='time') - ozone_gromos.mean(dim='time'))/o3_mls.mean(dim='time')
   # rel_diff_somora_mls = 100*(o3_mls.mean(dim='time') - ozone_somora.mean(dim='time'))/o3_mls.mean(dim='time')

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(9, 6))

    ozone_somora.where(mr_somora > 0.8, drop=False).mean(dim='time').plot(
        y='o3_p', ax=axs[0], color='red', label='SOMORA')
    ozone_gromos.where(mr_gromos > 0.8, drop=False).mean(dim='time').plot(
        y='o3_p', ax=axs[0], color='blue', label='GROMOS')
    o3_mls.mean(dim='time').plot(y='p', ax=axs[0], color='black', label='MLS')

    gromos_plus = ozone_gromos.mean(dim='time') + gromos_error.mean(dim='time')
    gromos_minus = ozone_gromos.mean(
        dim='time') - gromos_error.mean(dim='time')
    somora_plus = ozone_somora.mean(dim='time') + somora_error.mean(dim='time')
    somora_minus = ozone_somora.mean(
        dim='time') - somora_error.mean(dim='time')

    gromos_minus.plot(
        y='o3_p', ax=axs[0], ls='--', color='blue', alpha=0.4, label='mean obs error')
    gromos_plus.plot(y='o3_p', ax=axs[0], ls='--', color='blue', alpha=0.4)

    somora_plus.plot(y='o3_p', ax=axs[0], ls='--',
                     color='red', alpha=0.4, label='mean obs error')
    somora_minus.plot(y='o3_p', ax=axs[0], ls='--', color='red', alpha=0.4)

   # apriori_mean.mean(dim='time').plot(y='o3_p',label='GROMOS')
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].legend()
    axs[0].set_ylim(500, 1e-3)
    axs[0].set_ylabel('pressure [hPa]')
    axs[0].grid(axis='x', linewidth=0.5)

    rel_diff.plot(y='o3_p', ax=axs[1], color='green',
                  ls='--', alpha=0.4, label='(red-blue)/red')
    rel_diff.where(mr_mean_gromos > 0.8).plot(
        y='o3_p', ax=axs[1], color='green', label='(red-blue)/red')
  #  rel_diff_gromos_mls.plot(y='o3_p',ax=axs[1], color='blue', label='(MLS-GRO)/MLS')
  #  rel_diff_somora_mls.plot(y='o3_p',ax=axs[1], color='red', label='(MLS-SOM)/MLS')
    axs[1].set_xlim(-60, 60)
    axs[1].set_xlabel('relative difference [%]')
    axs[1].legend()
    axs[1].grid(axis='x', linewidth=0.5)

    for a in axs:
        a.grid(which='both', axis='y', linewidth=0.5)
    plt.suptitle('Ozone comparison with ' + str(len(date)) + ' days ' +
                 pd.to_datetime(ozone_somora.time.mean().data).strftime('%Y-%m-%d %H:%M'))
    plt.tight_layout()

    # save mean profiles for these dates
    o3_mls.mean(dim='time').to_netcdf(somora.level2_folder+'/' +
                                      'mls_mean_o3_'+date.mean().strftime('%Y-%m-%d')+'.nc')
    ozone_gromos.mean(dim='time').to_netcdf(
        somora.level2_folder+'/'+'gromos_mean_o3_'+date.mean().strftime('%Y-%m-%d')+'.nc')
    ozone_somora.mean(dim='time').to_netcdf(
        somora.level2_folder+'/'+'somora_mean_o3_'+date.mean().strftime('%Y-%m-%d')+'.nc')

    fig.savefig(somora.level2_folder+'/'+'ozone_comparison_' +
                pd.to_datetime(ozone_somora.time.mean().data).strftime('%Y-%m-%d')+'.pdf')

if plot_all:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + '_plot_all_test_polyfit2'
    GROSOM_library.plot_O3_all(level2_dataset, outname)

if plot_o3_ts:
    ozone = level2_dataset['AC240'].isel(o3_lat=0, o3_lon=0)
    o3 = ozone.o3_x
    mr = ozone.o3_mr.data
    o3['altitude'] = ozone.o3_z / 1e3
    # o3_test = xr.DataArray(
    #     data=ozone.o3_x.data,
    #     dims={'time','altitude'},
    #     coords=dict(
    #         time=ozone.time.values,
    #         altitude=(['time','altitude'],ozone.o3_z.values)
    #     )
    # )
   # o3 = o3.assign_coords({'altitude':ozone.o3_z})
    o3_hourly = o3.resample(time='8H', skipna=True).nearest(tolerance='1H')
    o3.coords['o3_p'] = o3.coords['o3_p']/100
    o3.data = o3.data*1e6
   # o3=o3.swap_dims({'o3_p':'altitude'})
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
  #  o3_hourly = o3_hourly.swap_dims({'o3_p':'geometric_height'})
   # o3_hourly['geometric_height'] = o3_hourly.o3_z
    #o3.plot(x='time', y='altitude')
    pl = o3.where(mr > 0.8).resample(time='4H', skipna=True).mean().plot(
        x='time',
        y='o3_p',
        vmin=0,
        vmax=9,
        cmap=colormap,
        yscale='log',
        # linewidth=0,
        # rasterized=True,
        cbar_kwargs={"label": "ozone [PPM]"}
    )
    pl.set_edgecolor('face')
    ax.set_ylim(0.01, 500)
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_yscale('log')

    ax.set_ylabel('Pressure [hPa]')
    plt.tight_layout()
    # o3.plot.imshow(x='time')
    fig.savefig(instrument.level2_folder+'/'+instrument.basename_plot_level2 +
                instrument.datestr+ex+'_ozone_ts_mr.pdf', dpi=500)
    if save_o3:
        o3_ds = ozone.get([
            'o3_x','o3_xa','o3_mr','o3_eo','o3_es','o3_avkm','o3_z','o3_fwhm', 'o3_offset',
            'median_noise','oem_diagnostics','obs_za','obs_aa','obs_lat','obs_lon','obs_alt']
            ) 
        o3_ds.to_netcdf(instrument.level2_folder+'/'+instrument.datestr+ex+'_ozone_ts_mr.nc')

if plot_o3_diff_waccm:
    # filename_waccm = '/storage/nas/MW/scratch/sauvageat/InputsRetrievals/waccm_o3_climatology.nc'
    # ds_waccm = xr.open_dataset(
    #         filename_waccm,
    #         decode_times=True,
    #         decode_coords=True,
    #         use_cftime=False,
    #     )

    # waccm_hourly = ds_waccm.isel(time=np.arange(0,90)).resample(time=1/24)
    
    ozone = level2_dataset['AC240'].isel(o3_lat=0, o3_lon=0)
    o3 = ozone.o3_x
    mr = ozone.o3_mr.data
    o3['altitude'] = ozone.o3_z / 1e3
    # o3_test = xr.DataArray(
    #     data=ozone.o3_x.data,
    #     dims={'time','altitude'},
    #     coords=dict(
    #         time=ozone.time.values,
    #         altitude=(['time','altitude'],ozone.o3_z.values)
    #     )
    # )
   # o3 = o3.assign_coords({'altitude':ozone.o3_z})
    o3_hourly = o3.resample(time='1H', skipna=True).nearest(tolerance='1H')

    o3_weekly= o3.resample(time='1d', skipna=True).nearest(tolerance='1H')

    o3.coords['o3_p'] = o3.coords['o3_p']/100
    o3.data = o3.data*1e6
   # o3=o3.swap_dims({'o3_p':'altitude'})
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
  #  o3_hourly = o3_hourly.swap_dims({'o3_p':'geometric_height'})
   # o3_hourly['geometric_height'] = o3_hourly.o3_z
    #o3.plot(x='time', y='altitude')
    pl = (o3_hourly).where(mr > 0.8).resample(time='4H', skipna=True).mean().plot(
        x='time',
        y='o3_p',
        vmin=0,
        vmax=9,
        cmap=colormap,
        yscale='log',
        # linewidth=0,
        # rasterized=True,
        cbar_kwargs={"label": "ozone [PPM]"}
    )
    pl.set_edgecolor('face')
    ax.set_ylim(0.01, 500)
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_yscale('log')

    ax.set_ylabel('Pressure [hPa]')
    plt.tight_layout()
    # o3.plot.imshow(x='time')
    fig.savefig(instrument.level2_folder+'/'+instrument.basename_plot_level2 +
                instrument.datestr+ex+'_ozone_ts_mr.pdf', dpi=500)

if plot_all_mopi5:
    for s in spectros:
        outname = basename_lvl2+'/'+s+'_'+instrument.datestr + '_plot_fascod'+ex
        mopi5_library.plot_O3_sel_paper(
            level2_dataset,
            outname,
            spectro=s,
            cycles=np.arange(0,15)
        )

if read_waccm_clim:
    filename_dt = os.path.join('/home/eric/Documents/PhD/GROSOM/InputsRetrievals', 'O3_apriori_daily_waccm_Bern_daylight.mat')
    filename_nt = os.path.join('/home/eric/Documents/PhD/GROSOM/InputsRetrievals', 'O3_apriori_daily_waccm2x_Bern_nighttime.mat')

    dt = scipy.io.loadmat(filename_dt)
    nt = scipy.io.loadmat(filename_nt)


    o3_dt = dt['O3_mean_smooth']
    o3_std_dt = dt['O3_std_smooth']
    p_dt = dt['p_mean_smooth']

    o3_nt = nt['O3_mean_smooth']
    o3_std_nt = nt['O3_std_smooth']
    p_nt = nt['p_mean_smooth']

    o3 = [dt['O3_mean_smooth'], nt['O3_mean_smooth']]
    o3_std = [dt['O3_std_smooth'], nt['O3_std_smooth']]
    p = [dt['p_mean_smooth'], nt['p_mean_smooth']]

    doy = np.arange(1,367)

    ds_o3_all = xr.Dataset(
        data_vars=dict(
            o3=(['tod','level', 'time'], o3),
            o3_std=(['tod','level', 'time'], o3_std),
            p=(['tod','level', 'time'], p),
        ),
        coords=dict(
            tod=['day', 'night'],
            time=doy,
            level=np.arange(0,66)
        ),
        attrs=dict(description='ozone climatology from waccm')
    )
    ds_o3_all.time.attrs['description'] = 'day of year'
    ds_o3_all.tod.attrs['description'] = 'time of day: daytime [0] or nighttime [1]'

    ds_o3_all.to_netcdf('waccm_o3_climatology.nc')

if compare_MLS:
    MLS_basename = '/home/esauvageat/Documents/AuraMLS/'
    filename_MLS = os.path.join(MLS_basename, 'aura-ozone-at-Bern.mat')
    mls_data = scipy.io.loadmat(filename_MLS)

    ozone_mls = mls_data['o3']
    p_mls = mls_data['p'][0]
    time_mls = mls_data['tm'][0]
    datetime_mls = []
    for i, t in enumerate(time_mls):
        datetime_mls.append(datetime.datetime.fromordinal(
            int(t)) + datetime.timedelta(days=time_mls[i] % 1) - datetime.timedelta(days=366))

    ds_mls = xr.Dataset(
        data_vars=dict(
            o3=(['time', 'p'], ozone_mls)
        ),
        coords=dict(
            lon=mls_data['longitude2'][0],
            lat=mls_data['latitude2'][0],
            time=datetime_mls,
            p=p_mls
        ),
        attrs=dict(description='ozone time series at bern')
    )
    #ds_mls = xr.decode_cf(ds_mls)
    #ds_mls.time.encoding['units'] = 'seconds since 1970-01-01 00:00:00'
    #ds_mls.time.encoding['calendar'] = "proleptic_gregorian"

    #ds_mls.to_netcdf('/home/esauvageat/Documents/AuraMLS/ozone_bern_ts.nc', format='NETCDF4')

    o3_mls = ds_mls.sel(time=slice("2019-01-30", "2019-06-22")).o3

    monthly_mls = o3_mls.resample(time='1m', skipna=True).mean()

    fig, ax = plt.subplots(1, 1)
    #monthly_mls.plot(x='time', y='p', ax=ax ,vmin=0, vmax=9).resample(time='24H', skipna=True).mean()
    pl = o3_mls.plot(
        x='time',
        y='p',
        ax=ax,
        vmin=0,
        vmax=10,
        yscale='log',
        linewidth=0,
        rasterized=True,
        cmap=colormap
    )
   # pl.set_edgecolor('face')
    # ax.set_yscale('log')
    ax.set_ylim(0.01, 500)
    ax.invert_yaxis()
    ax.set_ylabel('P [hPa]')
    plt.tight_layout()
    fig.savefig(instrument.level2_folder +
                '/ozone_mls_01-02-2019_mr'+'.pdf', dpi=500)


if plot_selected:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + '_plot_sel_polyfit2'
    GROSOM_library.plot_O3_all(
        level2_dataset,
        outname,
        spectro='AC240',
        cycles=[1,6,9]
    )

if plot_fshift:
    outname = instrument.level2_folder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + '_plot_fshift.pdf'
    fshift = level2_dataset['AC240'].isel(freq_shift_grid1=0).freq_shift_x

    fshift = fshift*1e-3
    fig = plt.figure(num=1)
    ax = fig.subplots(1)
    fshift.resample(time='4H', skipna=True).mean().plot(
        ax=ax
    )
    ax.set_ylabel('fshift [kHz]')
    ax.set_title('fshift '+instrument_name)
    plt.tight_layout()
    fig.savefig(outname)

if compare_MERRA2:
    merra2_basename = '/storage/tub/atmosphere/MERRA2/BRN/'
    filename_merra2 = [
        # os.path.join(merra2_basename,'MERRA2_BRN_2016_01_diagnostic.h5'),
        os.path.join(merra2_basename, 'MERRA2_BRN_2017_08_diagnostic.h5'),
        os.path.join(merra2_basename, 'MERRA2_BRN_2017_09_diagnostic.h5'),
        os.path.join(merra2_basename, 'MERRA2_BRN_2017_10_diagnostic.h5'),
        os.path.join(merra2_basename, 'MERRA2_BRN_2017_11_diagnostic.h5'),
        os.path.join(merra2_basename, 'MERRA2_BRN_2017_12_diagnostic.h5')
    ]

    o3_merra2_tot = xr.Dataset()
    counter = 0
    for f in filename_merra2:
        merra2_info = xr.open_dataset(
            f,
            group='info',
            decode_times=False,
            decode_coords=True,
            # use_cftime=True,
        )
        #merra2_info.time.attrs = attrs
        merra2_info['time'] = merra2_info.time - merra2_info.time[0]
    #   # construct time vector

        time_merra2 = []
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
            # use_cftime=True,
        )

        o3_merra2 = merra2.O3

        o3_merra2 = o3_merra2.swap_dims({'phony_dim_6': 'altitude'})
        o3_merra2['altitude'] = merra2_decoded.alt.isel(phony_dim_1=0).data

        o3_merra2 = o3_merra2.swap_dims({'phony_dim_5': 'datetime'})
        o3_merra2['datetime'] = merra2_decoded.datetime.data

        o3_merra2.data = o3_merra2.data*1e6
        if counter == 0:
            o3_merra2_tot = o3_merra2
        else:
            o3_merra2_tot = xr.concat(
                [o3_merra2_tot, o3_merra2], dim='datetime')
        counter = counter + 1

    o3_merra2_tot = o3_merra2_tot.sel(
        datetime=slice("2017-08-15", "2017-12-31"))

    fig, ax = plt.subplots(1, 1)
    o3_merra2_tot.plot(x='datetime', y='altitude', ax=ax, vmin=0, vmax=15, cmap=colormap)
    ax.set_ylim(5, 65)
    # o3_merra2.assign_coords({'altitude':merra2_info.alt.isel(phony_dim_1=0)})
    plt.tight_layout()
    fig.savefig(instrument.level2_folder+'/ozone_ts_17_merra2.pdf')
# %%

if compare_ECMWF:
    ECMWF_folder = '/scratch/ECMWF/'
    counter = 0
    for d in date:
        ECMWF_file = os.path.join(
            ECMWF_folder, 'ecmwf_oper_v2_PAYERNE_'+d.strftime('%Y%m%d')+'.nc')

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
    # o3.plot.imshow(x='time')
    fig2.savefig(instrument.level2_folder+'/'+'ozone_ts_16_ecmwf_payerne.pdf')

    from retrievals.data.ecmwf import ECMWFLocationFileStore, levels

    ecmwf_prefix = 'ecmwf_oper_v2_BERN_%Y%m%d.nc'
    t1 = date[0]
    t2 = date[2]
    ecmwf_store = ECMWFLocationFileStore(ECMWF_folder, ecmwf_prefix)
    ds_ecmwf = (
        ecmwf_store.select_time(t1, t2, combine='by_coords')
        .mean(dim='time')
        .swap_dims({"level": "pressure"})
    )

    ds_ecmwf = read_add_geopotential_altitude(ds_ecmwf)
