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
    #date = datetime.date(2019,2,21)
    #date = pd.date_range(start='2019-01-03', end='2019-01-05')
    #meanTb_chunks = [100, 110, 120, 130, 140, 160, 180]


    #date = pd.date_range(start='2019-01-30', end='2019-06-18')

    #date = pd.date_range(start='2019-01-30', end='2019-02-22')
    #meanTb_chunks = [80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 190]

    #date = pd.date_range(start='2019-05-01', end='2019-05-04')
    # No U5303

    #date = pd.date_range(start='2019-04-25', end='2019-04-27')
    #meanTb_chunks = [105, 110, 115, 120, 130, 140, 160, 180, 200]


    #date = pd.date_range(start='2019-06-11', end='2019-06-15')
    #meanTb_chunks = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220]


    #date = pd.date_range(start='2019-03-12', end='2019-03-12')
    date = pd.date_range(start='2019-03-12', end='2019-03-12')
    # options are: 'TOD', 'TOD_harmo', 'classic' 'meanTb_harmo', or 'meanTb'

    integration_strategy = 'meanTb_harmo'
    int_time = 1

    plot_ts_Tb_Tsys = False
    df_bins=200e3
    date1b = pd.to_datetime(date[-1])
    #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/"

    plot_comparison = False
    plot_bias = True

    # Define the parameters for integration
    #TOD = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]
    #interval = np.ones(len(TOD))
    TOD = np.arange(24)
    interval = 0.5*np.ones(len(TOD))
    #TOD = [3, 9, 15, 21]

    classic = np.arange(1,24)

# %%
        
    basename_lvl1 = "/scratch/MOPI5/Level1/"
    basename_lvl2 = "/scratch/MOPI5/Level2/"
    #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    #calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time, ['AC240','USRP-A'])
    
    integration = mc.MOPI5_LvL2(date1b, basename_lvl1, integration_strategy, integration_strategy, integration_time=int_time)
    
    # Plotting part
    integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(no_flag=True, meta_data=False)

    if integration_strategy == 'meanTb' or integration_strategy == 'meanTb_harmo':
        identifier_plot = integrated_data['AC240'].coords['chunks'].data[1:].tolist() + [300]
        #identifier_plot = meanTb_chunks  + [300]
        idx_all = np.arange(0,len(identifier_plot))
        dimension=['chunks','channel_idx']
    else:
        dimension=['time','channel_idx']
        identifier_plot = TOD
        idx_all = np.arange(0,len(TOD))
    
    if plot_comparison:
        integration.compare_spectra_mopi5(
            dim=dimension[0], 
            idx=idx_all, 
            #idx=[0,1,2,3],
            save_plot = True, 
            #identifier=TOD,
            identifier=identifier_plot,
            with_corr = True,
            corr_band=False
        )
        integration.plot_time_min_comp()

        integration.compare_spectra_binned_interp_mopi5(
                dim=dimension[0], 
                idx=idx_all, 
                #idx=[0,1,2,3],
                spectrometers=integration.spectrometers,
                save_plot = True, 
                #identifier=TOD,
                identifier=identifier_plot,
                use_basis='U5303'
            )
        integration.compare_spectra_binned_interp_mopi5(
                dim=dimension[0], 
                #idx=np.arange(0,len(meanTb_chunks)+1), 
                idx=idx_all,
                spectrometers=integration.spectrometers,
                save_plot = True, 
                identifier=identifier_plot,
                use_basis='U5303',
                corrected=True
            ) 
    
    if plot_bias:
        figures=list()
        fig = plt.figure()
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)
        # ax1 = fig.add_subplot(3,1,1)
        # ax2 = fig.add_subplot(3,1,2)
        # ax3 = fig.add_subplot(3,1,3)
        end_dates = [datetime.date(2019,1,5), datetime.date(2019,2,22), datetime.date(2019,4,27)]
        #end_dates = pd.date_range(start='2019-01-03', end='2019-04-30') 
        #end_dates = [datetime.date(2019,1,5), datetime.date(2019,2,22), datetime.date(2019,1,30) ]
        #end_dates = [datetime.date(2019,2,12), datetime.date(2019,2,13), datetime.date(2019,3,12), datetime.date(2019,6,13),datetime.date(2019,6,14)]
        monthly_color = ['magenta', 'blue', 'cyan', 'orange', 'red']
        month_name = ['Jan', 'Feb', 'Mar', 'Apr']
        # month_name = ['Jan', 'Feb', 'Mar', 'Apr']
        # month = [pd.date_range(start='2019-01-01', end='2019-01-31'),pd.date_range(start='2019-02-01', end='2019-02-28'),pd.date_range(start='2019-03-01', end='2019-03-31'),pd.date_range(start='2019-04-01', end='2019-04-30')]

        size=3
        #integrated_data_monthly = xr.Dataset()
        #for m in month:
        count = 0
        for d in end_dates:
            try:
                integration = mc.MOPI5_LvL2(d, basename_lvl1, integration_strategy, integration_strategy, integration_time=int_time)
                integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(no_flag=True, meta_data=False)
                color = monthly_color[d.month-1]
                #color = monthly_color[d.day-11]
                # ax1.scatter(integrated_data['AC240'].time_min, integrated_data['AC240'].mean_hot_counts, s=size)
                # ax2.scatter(integrated_data['U5303'].time_min, integrated_data['U5303'].mean_hot_counts, s=size)
                # ax3.scatter(integrated_data['USRP-A'].time_min, integrated_data['USRP-A'].mean_hot_counts, s=size)
                s='AC240'    
                scatter = ax1.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].bias_Tb_lc.data, c=color, s=size)
                ax1.set_ylabel(r'$\Delta T_b$ [K]')
                ax1.set_ylim(0,2)
                ax1.set_title('Line center bias')
                ax1.set_xlabel('Mean $T_b$ [K]')
                ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope*1e9, c=color, s=size)
                ax3.set_ylabel('m [K/GHz]') 
                ax3.set_ylim(-0.5,1)
                ax3.set_title('Slope difference')
                ax3.set_xlabel('Mean $T_b$ [K]')
                ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].bias_Tb_lc_corr.data, c=color, s=size)
                ax2.set_ylabel(r'$\Delta T_b$ [K]')
                ax2.set_xlabel('Mean $T_b$ [K]')
                ax2.set_title('$\Delta T_b$ after trop corr')
                ax2.set_ylim(0,2)
            except:
                print('no data for :', d)
                pass           
        
        # handles, labeling = scatter.legend_elements(
        #     prop='colors',
        # )
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='.', color='w', markerfacecolor='magenta', label=month_name[0], markersize=size+5),
            Line2D([0], [0], marker='.', color='w', markerfacecolor='blue', label=month_name[1], markersize=size+5),
            Line2D([0], [0], marker='.', color='w', markerfacecolor='cyan', label=month_name[2], markersize=size+5),
            Line2D([0], [0], marker='.', color='w', markerfacecolor='orange', label=month_name[3], markersize=size+5)
        ]
        legend_elements = [
            Line2D([0], [0], marker='.', color='w', markerfacecolor='magenta', label=month_name[0], markersize=size+5),
            Line2D([0], [0], marker='.', color='w', markerfacecolor='blue', label=month_name[1], markersize=size+5),
            Line2D([0], [0], marker='.', color='w', markerfacecolor='orange', label=month_name[3], markersize=size+5)
        ]

        #legend = ax1.legend(*scatter.legend_elements(prop='colors'), month_name, fontsize='xx-small',loc=1, title='Month')
        #ax1.add_artist(legend)
        ax1.legend(handles=legend_elements, fontsize='small',loc=1, title='Month')
        ax1.grid()
        ax2.grid()
        ax3.grid()
        fig.suptitle('Bias '+s+' vs U5303')
        #fig.suptitle('Mean hot counts')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        
        figures.append(fig)
        save_single_pdf(basename_lvl1+'bias_hourly_'+s+'_'+integration_strategy+'.pdf',figures)


# %%    