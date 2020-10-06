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

from retrievals import data

import mopi5_classes as mc

# %%

if __name__ == "__main__":
    instrument_name = "mopi5"
    #date = datetime.date(2019,2,1)
    date = pd.date_range(start='2019-01-03', end='2019-01-05')
    #date = pd.date_range(start='2019-01-30', end='2019-06-18')
    date = pd.date_range(start='2019-01-30', end='2019-02-22')

    #date = pd.date_range(start='2019-05-01', end='2019-05-04')
    date = pd.date_range(start='2019-04-25', end='2019-04-27')
    date = pd.date_range(start='2019-06-11', end='2019-06-15')
    #date = pd.date_range(start='2019-03-12', end='2019-03-12')
    # options are: 'TOD', 'TOD_harmo', 'classic' 'meanTb_harmo', or 'meanTb'
    integration_strategy = 'TOD_harmo'
    int_time = 3
    save_nc = True
    plot_ts_Tb_Tsys = False
    df_bins=200e3
    date1b = pd.to_datetime(date[-1])
    #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/"

    plot_comparison = False

    # Define the parameters for integration
    TOD = [3, 9, 15, 21]
    interval = [3, 3, 3, 3]
    #meanTb_chunks = [80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 190]
    meanTb_chunks = [100, 110, 120, 130, 140, 160, 180]
    meanTb_chunks = [105, 110, 115, 120, 140, 160, 180, 200]
    meanTb_chunks = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220]
    classic = np.arange(1,24)

    if integration_strategy == 'meanTb' or integration_strategy == 'meanTb_harmo':
        dimension=['chunks','channel_idx']
    else:
        dimension=['time','channel_idx']

# %%
        
    basename_lvl1 = "/scratch/MOPI5/Level1/"
    basename_lvl2 = "/scratch/MOPI5/Level2/"
    #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    #calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time, ['AC240','USRP-A'])
    
    #integration = mc.MOPI5_LvL2(date1b, basename_lvl1, integration_strategy, integration_strategy, integration_time=int_time)

    # Plotting part
    #integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(no_flag=True, meta_data=False)

    if plot_comparison:
        integration.compare_spectra_mopi5(
            dim=dimension[0], 
            #idx=np.arange(0,len(meanTb_chunks)+1), 
            idx=[0,1,2,3],
            save_plot = True, 
            identifier=TOD,
            #identifier=meanTb_chunks+[300],
            with_corr = True
        )
        integration.plot_time_min_comp()

        integration.compare_spectra_binned_interp_mopi5(
                dim=dimension[0], 
                #idx=np.arange(0,len(meanTb_chunks)+1), 
                idx=[0,1,2,3],
                spectrometers=integration.spectrometers,
                save_plot = True, 
                identifier=TOD,
                #identifier=meanTb_chunks+[300],
                use_basis='U5303'
            ) 

    fig = plt.figure()
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    end_dates = [datetime.date(2019,1,5), datetime.date(2019,2,22), datetime.date(2019,4,27), datetime.date(2019,3,12), datetime.date(2019,6,15)]
    end_dates = pd.date_range(start='2019-02-01', end='2019-06-18') 
    #end_dates = [datetime.date(2019,2,12), datetime.date(2019,2,13), datetime.date(2019,3,12), datetime.date(2019,6,13),datetime.date(2019,6,14)]

    for d in end_dates:
        try:
            integration = mc.MOPI5_LvL2(d, basename_lvl1, integration_strategy, integration_strategy, integration_time=int_time)
            integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(no_flag=True, meta_data=False)
            for s in ['AC240']:
                ax1.plot(integrated_data[s].mean_Tb, integrated_data[s].bias_Tb_lc,'o',label=d.strftime('%Y-%m'))
                ax1.set_ylabel(r'$\Delta T_b$ [K]')
                ax2.plot(integrated_data[s].mean_Tb, integrated_data[s].slope*1e9,'o',label=d.strftime('%Y-%m'))
                ax2.set_ylabel('slope [K/GHz]')
                ax3.plot(integrated_data[s].mean_Tb, integrated_data[s].f_deltaTb0/1e9,'o',label=d.strftime('%Y-%m'))
                ax3.set_ylabel('f [GHz]')
                ax3.legend(fontsize='xx-small',loc=1)
        except:
            print('no data for :', d)
            pass

        

# %%