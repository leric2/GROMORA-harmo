#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Jun 15 2022

@author: eric

Scripts for the concatenation of the daily level 2 files from the GROMORA project.
It also adds the retrieval_quality flags before saving it to yearly netCDF file. 

This scripts also provides basic plotting capabilities of the main diagnostics quantities.
"""

import sys, os
from os.path import dirname, abspath, join

sys.path.append(join(dirname(sys.path[0]),'pyretrievals'))
sys.path.append(join(dirname(sys.path[0]),'retrieval'))

import datetime

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv

import GROMORA_library
# import mopi5_library
from gromora_utils import save_single_pdf

#from cmcrameri import cm
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.sans-serif": ["Free sans"]})

#load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')
#load_dotenv('/opt/arts/.env.stockhorn-arts24')
#load_dotenv('/opt/anaconda/.env.birg-arts24')

instrument_name = "GROMOS"

date = pd.date_range(start=sys.argv[1], end=sys.argv[2])
retrieval_strategy = sys.argv[3]
#date = datetime.date(2016,1,2)
#date = pd.date_range(start='2023-01-01', end='2023-01-01') 
#[pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=7)), pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=6))]


#############################################################################
# Concatenation of the daily files to single netCDF file for the period provided 
concatenate_and_add_L2_flags = True
save_residuals=False
save_level2 = True

# Base diagnostic plot for the day provided.
plot_selected = False

# raw plot of o3 time series and other diagnostics quantities
plot_o3_ts = True
plot_fshift = True
plot_cost = True
plot_polyfit = True
plot_sinefit = True

add_opacity=False

integration_strategy = 'classic'
spectros = ['AC240']
spectro = spectros[0]
int_time = 1

if retrieval_strategy=='consolidated':
    ex = '_v3'
    plotfolder = '/storage/atmosphere/instruments/gromos/level2/GROMORA/v3/'
elif retrieval_strategy == 'oper':
    ex = '_oper'  
    plotfolder = '/storage/atmosphere/instruments/gromos/level2/GROMORA/oper/'
else:
    raise ValueError('Retrieval strategy not valid !')
    
new_L2 = True

# Deal with CET now ? Only for FB relevant
change2UTC = True

#plotfolder = '/scratch/GROSOM/Level2/GROMORA_retrievals_v2/'

cont_name = 'h2o_continuum_x' 

colormap = 'cividis'  # 'viridis' #, batlow_map cmap_crameri cividis

#############################################################################
#############################################################################
if instrument_name == "GROMOS":
    basename_lvl1 = "/storage/atmosphere/instruments/gromos/level1/GROMORA/"+str(date[0].year)
    #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    if new_L2:
        basename_lvl2 = "/storage/atmosphere/instruments/gromos/level2/GROMORA/v3/"+str(date[0].year)
    else:
        basename_lvl2 = "/storage/atmosphere/instruments/gromos/level2/GROMORA/v1/"+str(date[0].year)
    if spectro == 'AC240':
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
    import somora_classes as somora
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    if new_L2:
    #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
        basename_lvl2 = "/storage/atmosphere/instruments/somora/level2/v2/"+str(date[0].year)
    else:
        basename_lvl2 = "/storage/atmosphere/instruments/somora/level2/v1/"+str(date[0].year)
    instrument = somora.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )

level2_dataset = instrument.read_level2(
    spectrometers=spectros,
    extra_base=ex
)

F0 = instrument.observation_frequency

#############################################################################
#############################################################################
if concatenate_and_add_L2_flags:

    new_ds = level2_dataset[spectro] 

    #############################################################################
    # Adding the retrieval_quality flags to the concatenated level 2:
    good_data = xr.where((np.abs(new_ds.oem_diagnostics[:, 2] - 1) < instrument.cost_threshold(date.year[0], version=3) ) & (np.abs(new_ds.poly_fit_x[:,0].data)<instrument.polyfit_threshold), True, False)

    if (spectro == 'FB') and change2UTC:
        # FB measured in CET ! 
        new_ds['time'] = new_ds['time'] - pd.Timedelta(1, 'hour')
        instrument.timezone = 'Z'
    
    new_ds['retrieval_quality'] = ('time', good_data.data.astype(int)) # good_data*1 # ('time', good_data.data.astype(int))
    new_ds['retrieval_quality'].attrs['standard_name'] = 'retrieval_quality'
    new_ds['retrieval_quality'].attrs['long_name'] = 'quality flag retrieval'
    new_ds['retrieval_quality'].attrs['units'] = '1'
    new_ds['retrieval_quality'].attrs['description'] = 'Quality flag of the retrievals from cost and polyfit term'
    # Combined flags 
    # diff_polyfit_1 = new_ds.poly_fit_x[1] - new_ds.poly_fit_x[1].mean(dim='time')
    # diff_polyfit_2 = new_ds.poly_fit_x[2] - new_ds.poly_fit_x[2].mean(dim='time')

    #  good_data = xr.where(np.abs(diff_polyfit_2)<0.2, True, False)
    # new_ds.where(good_data).o3_x.isel(o3_p=15).plot()

    #############################################################################
    # Adding further attributes
    new_ds['o3_x'].attrs['valid_min'] = 0
    new_ds['o3_x'].attrs['valid_max'] = 50e-6
    
    # Encoding the time coordinates properly:
    new_ds.time.attrs['standard_name'] = 'time'
    new_ds.time.encoding['units'] = 'days since 2000-01-01 00:00:00'
    new_ds.time.encoding['calendar'] = 'proleptic_gregorian' #'standard'
    new_ds.time.attrs['timezone'] = instrument.timezone
    new_ds.time.attrs['description'] = 'mean time recorded at the beginning of all sky measurements during this integration cycle'
    
    #############################################################################
    # For saving the residuals, for instance for baseline analysis
    if save_residuals:
        residual = new_ds.get(['y', 'yf', 'bad_channels', 'y_baseline'])
        residual['res'] = residual.y - residual.yf
        residual.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr+ex+'_residuals.nc')
    
    # If not, we remove the observation vector from the concatenated level 2 (for space)
    ozone = new_ds.drop_dims(['f']) #drop_vars(['y', 'yf', 'bad_channels', 'y_baseline'])
    if save_level2:
        ozone.to_netcdf(plotfolder+'/'+instrument_name+'_'+instrument.datestr+'_'+spectro+ex+'.nc' , format='NETCDF4', unlimited_dims='time')

#############################################################################
#############################################################################
if plot_selected:
    outname = plotfolder+'/'+instrument.basename_plot_level2 + \
        instrument.datestr + ex + '_plot_sel_v2'
    

    instrument.plot_ozone_sel(
        level2_dataset,
        outname,
        spectro=spectros[0],
        cycles=[0,8,16],
        altitude = False,
        add_baselines = True, 
        to_ppm = 1e6  
    )

#############################################################################
#############################################################################
if plot_o3_ts:
    ozone = level2_dataset[spectro] 
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
    pl = o3.resample(time='1H', skipna=True).mean().plot(
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
    #axs.set_ylim(0.01, 500)
    axs.set_yticks([])
    axs.invert_yaxis()
    axs.set_yscale('log')

    axs.set_ylabel('Pressure [hPa]')
    plt.tight_layout()
    # o3.plot.imshow(x='time')
    fig2.savefig(plotfolder+'/'+instrument.basename_plot_level2 +
                instrument.datestr+ex+'_ozone_ts_mr.pdf', dpi=500)

#############################################################################
#############################################################################
if plot_polyfit:
    for s in spectros:
        fig, axes = plt.subplots(nrows=4, ncols=1,sharex=True, figsize=(15, 10))
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
        
        fig.tight_layout(rect=[0, 0.03, 1, 0.99])
        fig.savefig(plotfolder+'/'+instrument.basename_plot_level2 +
                        instrument.datestr+ex+'_polyfit.pdf', dpi=500)

#############################################################################
#############################################################################
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
        instrument.datestr+ex+'_sinefit.pdf', dpi=500)      

#############################################################################
#############################################################################
if plot_cost:
    for s in spectros:
        fig, axsco = plt.subplots(nrows=2, ncols=1,sharex=True, figsize=(15, 10))
        end_cost = level2_dataset[s].oem_diagnostics[:, 3]
        noise_input = level2_dataset[s].median_noise
        end_cost.plot(marker='.', ax=axsco[0])
        noise_input.plot(marker='.',ax=axsco[1])

          #  end_cost.plot(ax=axs, ylim=(0.75,8))
        fig.savefig(plotfolder+'/'+instrument.basename_plot_level2 + instrument.datestr+ex+'_end_cost.pdf', dpi=500)

