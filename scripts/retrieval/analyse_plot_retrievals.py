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

import sys, os
from os.path import dirname, abspath, join

sys.path.append(join(dirname(sys.path[0]),'pyretrievals'))
sys.path.append(join(dirname(sys.path[0]),'retrieval'))

import datetime
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
from xarray.backends import file_manager

import GROMORA_library
# import mopi5_library
from gromora_utils import save_single_pdf

import matplotlib
matplotlib.use('pdf')

#from cmcrameri import cm
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.sans-serif": ["Free sans"]})

# plt.rcParams['xtick.labelsize'] = 24
# plt.rcParams['ytick.labelsize'] = 24
# plt.rcParams['axes.titlesize'] = 24
#load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')
#load_dotenv('/opt/anaconda/.env.birg-arts24')
#load_dotenv('/opt/arts/.env.stockhorn-arts24')
# ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
# ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
# ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']
# from apriori_data_GROSOM import read_add_geopotential_altitude
# if __name__ == "__main__":

instrument_name = "GROMOS"

# date = pd.date_range(start='2019-01-03', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

#date = pd.date_range(start=sys.argv[1], end=sys.argv[2])
#date = pd.date_range(start='2011-01-01', end='2011-12-31')
#date = datetime.date(2016,1,2)
date = [pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=4)),pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=3))]
# pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=3))]

int_time = 1

plot_cycle = False
df_bins = 200e3

plot_all = False
plot_all_mopi5 = False
plot_o3_ts = False
save_o3= False
save_o3_only = False
save_residuals=False
plot_selected = True
plot_selected_nicer = False
plot_fshift = False
save_fshift = False
plot_cost = False
plot_o3_diff_waccm = False
read_waccm_clim = False
compare = False
plot_polyfit = False
plot_sinefit = False
plot_MLS = False
add_L2_flags = False
compare_spectra = False
add_opacity=False

extract_fgrid = False

#filename_opacity = '/scratch/GROSOM/Level2/opacities/'+instrument_name+'/'+instrument_name+'_opacity_'+str(date[0].year)+'.nc'
#  filename_opacity_new = '/scratch/GROSOM/Level1/'+instrument_name+'_level1b_v2_all.nc'

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

ex=''

ex = '_waccm_monthly_continuum'
ex = '_waccm_cov_yearly_sza' #
ex = '_waccm_continuum'
ex = '_waccm_monthly_scaled_h2o'
ex = '_gromosAP_scaled_h2o'
ex = '_waccm_low_alt_dx10_nonWinCorr'
ex = '_gromosAP_low_alt'

ex = '_sinefit_optimized'
ex = '_waccm_low_alt_dx10'
ex = '_rect_SB'
ex = '_oper'
# ex = '_waccm_low_alt'

new_L2 = True

if new_L2:
    plotfolder = '/scratch/GROSOM/Level2/GROMORA_retrievals_v2/'
    plotfolder = '/storage/tub/instruments/gromos/level2/GROMORA/oper/'
    cont_name = 'h2o_continuum_x' 
else:
    cont_name = 'h2o_pwr98_x'
    plotfolder = '/scratch/GROSOM/Level2/GROMORA_waccm/'
# %%

colormap = 'cividis'  # 'viridis' #, batlow_map cmap_crameri cividis

def read_opacity(filename, date):
    ds_opacity = xr.open_dataset(
        filename,
        group='spectrometer1',
        decode_times=True,
        decode_coords=True,
        # use_cftime=True,
    )

    ds_opacity = ds_opacity.sel(time=slice(date[0].strftime('%Y-%m-%d'), date[-1].strftime('%Y-%m-%d')))
    return ds_opacity

def read_level1b(filename):
    ds = xr.open_dataset(
        filename,
        decode_times=True,
        decode_coords=True,
        # use_cftime=True,
    )
    return ds

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
    basename_lvl1 = "/storage/tub/instruments/gromos/level1/GROMORA/"+str(date[0].year)
    #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    if new_L2:
        basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v2/"+str(date[0].year)
    else:
        basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v1/"+str(date[0].year)
    if spectros[0] == 'AC240':
        import gromos_classes as gromos_cl
        instrument = gromos_cl.GROMOS_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
        )
    else:
        import gromos_FB_classes as gromos_cl
        instrument = gromos_cl.GROMOS_FB_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
elif instrument_name == "SOMORA":
    import somora_classes as sm
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    if new_L2:
    #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
        basename_lvl2 = "/storage/tub/instruments/somora/level2/v2/"+str(date[0].year)
    else:
        basename_lvl2 = "/storage/tub/instruments/somora/level2/v1/"+str(date[0].year)
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
    import gromos_FB_classes as gromos_cl_FB
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    #basename_lvl2 = "/storage/tub/instruments/somora/level2/v1/"+str(date[0].year)
    basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v2/"+str(date[0].year)
    somora = gromos_cl_FB.GROMOS_FB_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
    import gromos_classes as gromos_cl
    basename_lvl1 = "/scratch/GROSOM/Level1/GROMOS/"
    basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v2/"+str(date[0].year)
    import gromos_classes as gromos_cl
    gromos = gromos_cl.GROMOS_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
        )

if instrument_name == "compare":
    ex='_waccm_low_alt'
    level2_FB = somora.read_level2(
        spectrometers=['FB'],
        extra_base='_v2'
    )
    level2_gromos = gromos.read_level2(
        spectrometers=['AC240'],
        extra_base='_v2'
    )
    F0 = somora.observation_frequency

else:
    # Plotting part
    level2_dataset = instrument.read_level2(
        spectrometers=spectros,
        extra_base=ex
    )
    # level2_dataset = instrument.read_level2_multidays(
    #     spectrometers=spectros,
    #     extra_base=ex
    # )
    #print(level2_dataset)
    F0 = instrument.observation_frequency

if compare_spectra:
    fig, axs = plt.subplots(nrows=1, ncols=1,sharex=True, figsize=(15, 10))

    level2_gromos['AC240'].y.mean(dim='time').plot(ax = axs)
    level2_FB['FB'].y.mean(dim='time').plot(ax = axs)
    plt.show()
    fig.savefig(plotfolder+'/'+gromos.basename_plot_level2 +gromos.datestr+ex+'_daily_spectra_FB_FFT.pdf', dpi=500)

if extract_fgrid:
    date = datetime.date(2017, 10, 12)
    import somora_classes as sm
    basename_lvl1 = "/storage/tub/instruments/somora/level1/v1/"+str(date.year)
    #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    basename_lvl2 = "/storage/tub/instruments/somora/level2/v1/"+str(date.year)
    instrument = sm.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
    integrated_data, integrated_flags, integrated_meteo = instrument.read_level1b(
    no_flag=True, meta_data=False)
    integrated_data['AC240'].isel(time=0).frequencies.to_netcdf('/scratch/GROSOM/Level1/frequency_grid.nc')
if plot_cost:
    if add_opacity:
        ds_opacity = read_opacity(filename_opacity, date)

    if instrument_name == 'mopi5':
        for s in spectros:
            fig, axsco = plt.subplots(nrows=1, ncols=1, sharex=True,figsize=(15, 10))
            end_cost = level2_dataset[s].oem_diagnostics[:, 3]
            axsco.plot(end_cost)
            fig.savefig(plotfolder+'/' +
                        instrument.datestr+s+ex+'_end_cost.pdf', dpi=500)
    else:
        for s in spectros:
            fig, axsco = plt.subplots(nrows=3, ncols=1,sharex=True, figsize=(15, 10))
            end_cost = level2_dataset[s].oem_diagnostics[:, 3]
            noise_input = level2_dataset[s].median_noise
            end_cost.plot(marker='.', ax=axsco[0])
            noise_input.plot(marker='.',ax=axsco[1])
            if add_opacity:
                ds_opacity.tropospheric_opacity.plot(marker='.',ax=axsco[2])
            #axs[1].set_ylim(0,1)
            axsco[0].set_xticks([])
            axsco[1].set_xticks([])
          #  end_cost.plot(ax=axs, ylim=(0.75,8))
            fig.savefig(plotfolder+'/'+instrument.basename_plot_level2 +
                        instrument.datestr+ex+'_end_cost.pdf', dpi=500)


if plot_polyfit:
    if add_opacity:
        ds_opacity = read_opacity(filename_opacity, date)

    if instrument_name=='compare':
        s = 'AC240'
        fig, axes = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(15, 10))
        polyfit_1 = level2_somora[s].poly_fit_x
        polyfit_2 = level2_gromos[s].poly_fit_x

        polyfit_1[0].plot(marker='.', ax=axes[0])
        polyfit_2[0].plot(marker='.', ax=axes[0])
        polyfit_1[1].plot(marker='.',ax=axes[1])
        polyfit_2[1].plot(marker='.', ax=axes[1])
        polyfit_1[2].plot(marker='.',ax=axes[2])
        polyfit_2[2].plot(marker='.', ax=axes[2])
        level2_somora[s][cont_name].plot(ax=axes[3])
        level2_gromos[s][cont_name].plot(ax=axes[3])
        axes[0].set_ylabel('polyfit')
        axes[1].set_ylabel('polyfit')
        axes[2].set_ylabel('polyfit')
        axes[3].set_ylabel(r'H$_2$O PWR98')
        axes[3].set_title('Continuum (retrieved)')
        for i in [0,1,2,3]:
            axes[i].set_xlabel('')
            axes[i].set_xticks([])
        
        if add_opacity:
            ds_opacity.tropospheric_opacity.plot(marker='.',ax=axes[4], color='k')
            axes[4].set_ylabel('opacity [-]')
            axes[4].set_xlabel('')
            axes[4].set_title('Tropospheric opacity (Ingold)')
      #  perc=100*(polyfit_1[0]-polyfit_2[0])/polyfit_1[0]
       # perc.plot(marker='.', ax=axes[3])
        axes[0].legend(('SOMORA', 'GROMOS')) 
        axes[0].set_xlabel('')
        axes[1].set_xlabel('')
        fig.savefig(plotfolder+'/'+gromos.basename_plot_level2 +
                        gromos.datestr+ex+'polyfit_comparison.pdf', dpi=500)
    else:
        for s in spectros:
            fig, axes = plt.subplots(nrows=5, ncols=1,sharex=True, figsize=(15, 10))
            polyfit = level2_dataset[s].poly_fit_x
            if new_L2:
                polyfit[:,0].plot(marker='.', ax=axes[0], color='k')
                polyfit[:,1].plot(marker='.',ax=axes[1], color='k')
                polyfit[:,2].plot(marker='.',ax=axes[2], color='k')
            else:
                polyfit[0].plot(marker='.', ax=axes[0], color='k')
                polyfit[1].plot(marker='.',ax=axes[1], color='k')
                polyfit[2].plot(marker='.',ax=axes[2], color='k')
            level2_dataset[s][cont_name].plot(ax=axes[3], color='k')
            axes[0].set_ylabel('polyfit')
            axes[1].set_ylabel('polyfit')
            axes[2].set_ylabel('polyfit')
            axes[3].set_ylabel(r'H$_2$O PWR98')
            axes[3].set_title('Continuum (retrieved)')

            # axes[0].set_ylim((-0.01,0.15))
            for i in [0,1,2,3]:
                axes[i].set_xlabel('')
                axes[i].set_xticks([])
            if add_opacity:
                ds_opacity.tropospheric_opacity.plot(marker='.',ax=axes[4], color='k')
            #ds_opacity.tropospheric_opacity_tc.plot(marker='.',ax=axes[3])
            #axes[3].set_ylim(0,20)
            #axs[1].set_ylim(0,1)
            #  end_cost.plot(ax=axs, ylim=(0.75,8))
            axes[4].set_ylabel('opacity [-]')
            axes[4].set_xlabel('')
            axes[4].set_title('Tropospheric opacity (Ingold)')
            fig.tight_layout(rect=[0, 0.03, 1, 0.99])
            fig.savefig(plotfolder+'/'+instrument.basename_plot_level2 +
                        instrument.datestr+ex+'polyfit.pdf', dpi=500)

if plot_sinefit:
    s = 'AC240'
     #level2_dataset[s]['time'] = pd.to_datetime(level2_dataset[s]['time'].data)
    fig, axes = plt.subplots(nrows=3, ncols=1,sharex=True, figsize=(15, 10))
    level2_dataset[s].sine_fit_0_x[:,0].plot(marker='.', ax=axes[0], color='k')
    level2_dataset[s].sine_fit_0_x[:,1].plot(marker='.', ax=axes[0], color='r')
    level2_dataset[s].sine_fit_1_x[:,0].plot(marker='.', ax=axes[1], color='k')
    level2_dataset[s].sine_fit_1_x[:,1].plot(marker='.', ax=axes[1], color='r')
    level2_dataset[s].sine_fit_2_x[:,0].plot(marker='.', ax=axes[2], color='k')
    level2_dataset[s].sine_fit_2_x[:,1].plot(marker='.', ax=axes[2], color='r')

    axes[0].set_ylabel('sinefit')
    axes[1].set_ylabel('sinefit')
    axes[2].set_ylabel('sinefit')

    axes[0].set_title('Period ='+str(level2_dataset[s].sine_fit_0_x.attrs['period MHz'])+' MHz')
    axes[1].set_title('Period ='+str(level2_dataset[s].sine_fit_1_x.attrs['period MHz'])+' MHz')
    axes[2].set_title('Period ='+str(level2_dataset[s].sine_fit_2_x.attrs['period MHz'])+' MHz')

    for i in [0,1]:
        axes[i].set_xlabel('')
        #axes[i].set_xticks([])
    #ds_opacity.tropospheric_opacity_tc.plot(marker='.',ax=axes[3])
    #axes[3].set_ylim(0,20)
    #axs[1].set_ylim(0,1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.99])
    fig.savefig(plotfolder+'/'+instrument.basename_plot_level2 +
        instrument.datestr+ex+'sinefit.pdf', dpi=500)      

def avk_smooth_mls(ds, o3_mls_mean):
    avkm=ds.o3_avkm.mean(dim='time').values
    o3_p = ds.o3_p.values
    o3_p_mls = o3_mls_mean.p.values*100
    idx = np.argsort(o3_p_mls, kind='heapsort')
    o3_p_mls_sorted = o3_p_mls[idx]  
    o3_mls_sorted = o3_mls_mean[idx]

    interpolated_mls = np.interp(np.log(o3_p), np.log(o3_p_mls_sorted),o3_mls_sorted.values*1e-6)
    o3_convolved = ds.o3_xa.mean(dim='time').values  + np.matmul(avkm,(interpolated_mls - ds.o3_xa.mean(dim='time').values))
    return o3_convolved

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
    o3_mls_mean = o3_mls.mean(dim='time')
    ds = level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0)
    
    o3_mls_convolved = avk_smooth_mls(ds, o3_mls_mean)

    rel_diff = 100*(ozone_somora.mean(dim='time') -
                    ozone_gromos.mean(dim='time'))/ozone_somora.mean(dim='time')
   # rel_diff_gromos_mls = 100*(o3_mls.mean(dim='time') - ozone_gromos.mean(dim='time'))/o3_mls.mean(dim='time')
   # rel_diff_somora_mls = 100*(o3_mls.mean(dim='time') - ozone_somora.mean(dim='time'))/o3_mls.mean(dim='time')

    p_somora_mr = 1e-2*level2_somora['AC240'].isel(o3_lat=0, o3_lon=0).o3_p.data[np.mean(mr_somora,0)>=0.8]
    p_gromos_mr = 1e-2*level2_gromos['AC240'].isel(o3_lat=0, o3_lon=0).o3_p.data[np.mean(mr_gromos,0)>=0.8]

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(9, 6))

    # ozone_somora.where(mr_somora > 0.8).mean(dim='time').plot(
    #     y='o3_p', ax=axs[0], color='red', label='waccm')
    # ozone_gromos.where(mr_gromos > 0.8).mean(dim='time').plot(
    #     y='o3_p', ax=axs[0], color='blue', label='gromosAP')
    ozone_somora.mean(dim='time').plot(
        y='o3_p', ax=axs[0], color='red', ls='-', label='SOMORA')
    ozone_gromos.mean(dim='time').plot(
        y='o3_p', ax=axs[0], color='blue', ls='-', label='GROMOS')

    if plot_MLS:
        o3_mls.mean(dim='time').plot(y='p', ax=axs[0], color='black', label='MLS')

        axs[0].plot(o3_mls_convolved*1e6, ozone_gromos.o3_p.values, '--', color='black', label='MLS, convolved')

    axs[0].axhline(y=p_somora_mr[0],ls=':' ,color='r', lw=1)
    axs[0].axhline(y=p_somora_mr[-1],ls=':', color='r', lw=1)
    axs[0].axhline(y=p_gromos_mr[0],ls=':', color='b', lw=1)
    axs[0].axhline(y=p_gromos_mr[-1],ls=':', color='b', lw=1)

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
    #axs[1].set_xlim(-25, 25)
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

    fig.savefig(plotfolder+'/'+'ozone_comparison_'+pd.to_datetime(ozone_gromos.time.mean().data).strftime('%Y-%m-%d')+ex+'.pdf')

if plot_all:
    outname = plotfolder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + '_plot_all_test_polyfit2'
    GROMORA_library.plot_O3_all(level2_dataset, outname)

if add_L2_flags:
    new_ds = level2_dataset['AC240'] 
    # Absolute flags 
    good_data = xr.where((np.abs(new_ds.oem_diagnostics[:, 2] - 1) < instrument.cost_threshold(date.year[0]) ) & (np.abs(new_ds.poly_fit_x[:,0].data)<instrument.polyfit_threshold), True, False)
   # good_data = xr.where( np.abs(new_ds.oem_diagnostics[:, 3]-1)<0.1, True, False)
   # good_data = xr.where((np.abs(new_ds.oem_diagnostics[:, 3] - 1) <0.1 ), True, False) 
    
   # np.sum(np.abs(1-gromos.oem_diagnostics[:,2])<0.05)/(len(gromos.time))
  #  np.sum(np.abs(1-somora.oem_diagnostics[:,2])<0.05)/(len(somora.time))

    # Combined flags 
    # diff_polyfit_1 = new_ds.poly_fit_x[1] - new_ds.poly_fit_x[1].mean(dim='time')
    # diff_polyfit_2 = new_ds.poly_fit_x[2] - new_ds.poly_fit_x[2].mean(dim='time')

    #  good_data = xr.where(np.abs(diff_polyfit_2)<0.2, True, False)
    # new_ds.where(good_data).o3_x.isel(o3_p=15).plot()
    new_ds['o3_x'].attrs['valid_min'] = 0
    new_ds['o3_x'].attrs['valid_max'] = 50e-6

    new_ds['retrieval_quality'] = ('time', good_data.data.astype(int)) # good_data*1 # ('time', good_data.data.astype(int))
    new_ds['retrieval_quality'].attrs['standard_name'] = 'retrieval_quality'
    new_ds['retrieval_quality'].attrs['long_name'] = 'quality flag retrieval'
    new_ds['retrieval_quality'].attrs['units'] = '1'
    new_ds['retrieval_quality'].attrs['description'] = 'Quality flag of the retrievals from cost and polyfit term'
    
    new_ds['o3_x'].attrs['valid_min'] = 0
    new_ds['o3_x'].attrs['valid_max'] = 50e-6

    # new_ds['time'] = pd.to_datetime(level2_dataset['AC240'].time)
    new_ds.time.attrs['standard_name'] = 'time'
    new_ds.time.encoding['units'] = 'days since 2000-01-01 00:00:00'
    new_ds.time.encoding['calendar'] = 'standard'
    new_ds.time.attrs['timezone'] = 'Z'
    new_ds.time.attrs['description'] = 'mean time recorded at the beginning of all sky measurements during this integration cycle'

    if save_residuals:
        residual = new_ds.get(['y', 'yf', 'bad_channels', 'y_baseline'])
        residual['res'] = residual.y - residual.yf
        residual.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr+ex+'_residuals.nc')
    
    ozone = new_ds.drop_dims(['f']) #drop_vars(['y', 'yf', 'bad_channels', 'y_baseline'])
    ozone.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr+ex+'.nc' , format='NETCDF4', unlimited_dims='time')


if plot_o3_ts:
    ozone = level2_dataset['AC240'] 
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
    fig2, axs = plt.subplots(nrows=1, ncols=1, figsize=(9, 6))

  #  o3_hourly = o3_hourly.swap_dims({'o3_p':'geometric_height'})
   # o3_hourly['geometric_height'] = o3_hourly.o3_z
    #o3.plot(x='time', y='altitude')
    pl = o3.where(mr > 0.8).resample(time='1H', skipna=True).mean().plot(
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
    axs.set_ylim(0.01, 500)
    axs.set_yticks([])
    axs.invert_yaxis()
    axs.set_yscale('log')

    axs.set_ylabel('Pressure [hPa]')
    plt.tight_layout()
    # o3.plot.imshow(x='time')
    fig2.savefig(plotfolder+'/'+instrument.basename_plot_level2 +
                instrument.datestr+ex+'_ozone_ts_mr.pdf', dpi=500)

    if save_o3_only:
        o3_ds = ozone.get([
            'o3_x','o3_xa','o3_mr','o3_eo','o3_es','o3_avkm','o3_z','o3_fwhm', 'o3_offset',
            'median_noise','oem_diagnostics','obs_za','obs_aa','obs_lat','obs_lon','obs_alt']
            ) 
        o3_ds.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr+ex+'_ozone_ony.nc')
    if save_o3:
        o3_ds = ozone.get([
            'o3_x','o3_xa','o3_mr','o3_eo','o3_es','o3_avkm','o3_z','o3_fwhm', 'o3_offset','h2o_pwr98_x','freq_shift_x',
            'poly_fit_x','obs_time',
            'median_noise','oem_diagnostics','obs_za','obs_aa','obs_lat','obs_lon','obs_alt']
            ) 
        o3_ds.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr+ex+'_ozone.nc')
    if save_residuals:
        # residual = ozone.get([
        #     'y',
        #     'yf','y_baseline',
        #     'median_noise','oem_diagnostics','obs_za','obs_aa','obs_lat','obs_lon','obs_alt']
        #     ) 
        # residual  = np.convolve(ozone.y - ozone.yf, np.ones(4) / 4, mode="same")
        residual  = ozone.y - ozone.yf

        residual.rename('residuals')
        residual.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr+ex+'_residuals.nc')

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
    outname = plotfolder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + ex + '_plot_sel_polyfit2'
    
    if new_L2:
        instrument.plot_ozone_sel(
            level2_dataset,
            outname,
            spectro=spectros[0],
            cycles=[0,8,16],
            altitude = False,
            add_baselines = True, 
            to_ppm = 1e6  
        )
    else:
        GROMORA_library.plot_O3_all(
            level2_dataset,
            outname,
            spectro=spectros[0],
            cycles=[1,7,13,21]
        )

if plot_selected_nicer:
    outname = plotfolder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + ex + '_plot_sel_polyfit2'
    GROMORA_library.plot_O3_sel_nicer(
        level2_dataset,
        outname,
        spectro=spectros[0],
        cycles=[1,14,18]#cycles=[1,7,10,13,17,21]
    )


if plot_fshift:
    outname = plotfolder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + ex + '_plot_fshift.pdf'
    fshift = level2_dataset['AC240'].isel(f_shift_grid=0).freq_shift_x

    fshift = fshift*1e-3
    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(9, 6))

    fshift.plot(
        ax=ax
    )
    ax.set_ylabel('fshift [kHz]')
    #ax.set_ylim(-50,350)
    ax.set_title('fshift '+instrument_name)
    plt.tight_layout()
    fig.savefig(outname)

    if save_fshift:
        fshift.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr[0:4]+'_fshift.nc')

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
    ax.set_ylim(5, 75)
    # o3_merra2.assign_coords({'altitude':merra2_info.alt.isel(phony_dim_1=0)})
    plt.tight_layout()
    fig.savefig(instrument.level2_folder+'/ozone_ts_17_merra2.pdf')

if compare_ECMWF:
    Mair = 28.9644
    Mozone= 47.9982

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

    o3_ecmwf = ecmwf_ts.isel(loc=0).ozone_mass_mixing_ratio*Mair/Mozone
    o3_ecmwf['pressure'] = ecmwf_ts['pressure'].isel(loc=0)/100
    o3_ecmwf.data = o3_ecmwf.data * 1e6

    fig2 = plt.figure(num=1)
    ax = fig2.subplots(1)
    #o3.plot(x='time', y='altitude')
    o3_ecmwf.plot(
        x='time',
        y='pressure',
        vmin=0,
        vmax=9,
        cmap=colormap,
        cbar_kwargs={"label": "ozone [PPM]"}
    )
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylabel('P [hPa]')
    plt.tight_layout()
    # o3.plot.imshow(x='time')
    fig2.savefig(instrument.level2_folder+'/'+'ozone_ts_16_ecmwf_payerne.pdf')

    # from retrievals.data.ecmwf import ECMWFLocationFileStore, levels

    # ecmwf_prefix = 'ecmwf_oper_v2_BERN_%Y%m%d.nc'
    # t1 = date[0]
    # t2 = date[2]
    # ecmwf_store = ECMWFLocationFileStore(ECMWF_folder, ecmwf_prefix)
    # ds_ecmwf = (
    #     ecmwf_store.select_time(t1, t2, combine='by_coords')
    #     .mean(dim='time')
    #     .swap_dims({"level": "pressure"})
    # )

    # ds_ecmwf = read_add_geopotential_altitude(ds_ecmwf)
