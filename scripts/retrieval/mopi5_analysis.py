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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
# %%

if __name__ == "__main__":
    instrument_name = "mopi5"
    #date = datetime.date(2019,2,21)
    #date = pd.date_range(start='2019-01-03', end='2019-01-05')
    #meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]

    #date = pd.date_range(start='2019-01-30', end='2019-06-18')

    date = pd.date_range(start='2019-01-30', end='2019-02-22')
    meanTb_chunks = [80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 190]

    #date = pd.date_range(start='2019-05-01', end='2019-05-04')
    # No U5303

    #date = pd.date_range(start='2019-04-25', end='2019-04-27')
    #meanTb_chunks = [105, 110, 115, 120, 130, 160, 180, 200]


    #date = pd.date_range(start='2019-06-11', end='2019-06-15')
    #meanTb_chunks = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220]


    #date = pd.date_range(start='2019-03-12', end='2019-03-12')
    #date = pd.date_range(start='2019-02-22', end='2019-02-22')
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
    plot_fancy = True
    plot_bias = False
    plot_o3 = False

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
    
    integration = mc.MOPI5_LvL2(date1b, basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)
    
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

    if plot_fancy:
            integration.compare_spectra_mopi5(
            dim=dimension[0], 
            idx=[3], 
            #idx=[0,1,2,3],
            save_plot = True, 
            #identifier=TOD,
            with_corr = False,
            corr_band=False,
            title='Integrated spectra with $T_{B,mean}$ between 90 and 100K',
        )

            integration.compare_spectra_binned_interp_mopi5(
                dim=dimension[0],
                idx=idx_all, 
                spectrometers=['AC240','USRP-A'],
                save_plot = True, 
                use_basis='U5303',
                #identifier=TOD,
                identifier=identifier_plot,
                clean=True
        )

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
        #end_dates = [ datetime.date(2019,4,27)]
        #end_dates = [datetime.date(2019,2,12), datetime.date(2019,2,13), datetime.date(2019,3,12), datetime.date(2019,6,13),datetime.date(2019,6,14)]
        monthly_color = ['magenta', 'blue', 'cyan', 'orange', 'red']
        month_name = ['Jan', 'Feb', 'Mar', 'Apr','May','Jun']
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
                s='USRP-A'    
                scatter = ax1.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].bias_Tb_lc.data, c=color, s=size)
                ax1.set_ylabel(r'$\Delta T_b$ [K]')
                #ax1.set_ylim(-2,0)
                ax1.set_ylim(-1,1)
                ax1.set_title('Line center bias')
                ax1.set_xlabel('Mean $T_b$ [K]')
                ax3.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].slope*1e9, c=color, s=size)
                ax3.set_ylabel('m [K/GHz]') 
                ax3.set_ylim(-1,0.5)
                
                ax3.set_title('Slope difference')
                ax3.set_xlabel('Mean $T_b$ [K]')
                ax2.scatter(integrated_data[s].mean_Tb.data, integrated_data[s].bias_Tb_lc_corr.data, c=color, s=size)
                ax2.set_ylabel(r'$\Delta T_b$ [K]')
                ax2.set_xlabel('Mean $T_b$ [K]')
                ax2.set_title('$\Delta T_b$ after trop corr')
                #ax2.set_ylim(-2,0)
                ax2.set_ylim(-1,1)
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
        ax1.legend(handles=legend_elements, fontsize='small',loc=4, title='Month')
        ax1.grid()
        ax2.grid()
        ax3.grid()
        fig.suptitle('Bias '+s+' vs U5303')
        #fig.suptitle('Mean hot counts')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        
        figures.append(fig)
        save_single_pdf(basename_lvl1+'bias_binned_'+s+'_'+integration_strategy+'.pdf',figures)

    if plot_o3:

        figures2=list()
        #fig, axs = plt.subplots(nrows=3, ncols=5, sharex=True, sharey=True)
        fig = plt.figure()
        fig.subplots_adjust(hspace=0.4, wspace=0.4)
        
        fig_mr = plt.figure()
        fig_mr.subplots_adjust(hspace=0.4, wspace=0.4)

        spectro_lvl2 = integration.spectrometers

        level2_data = integration.read_level2(spectrometers=spectro_lvl2, extra_base='_all')

        color_spectro = {'AC240':'tab:orange', 'USRP-A':'tab:green', 'U5303':'tab:blue'}

        for spectro in spectro_lvl2:
            for i in range(len(level2_data[spectro].observation)):
                ax = fig.add_subplot(3,5,i+1)
                ax2 = fig_mr.add_subplot(3,5,i+1)

                o3 = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_x
                o3_apriori = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_xa
                o3_z = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_z
                o3_p = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_p
                mr = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_mr
                #error = lvl2[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_eo +  lvl2[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_es
                error = np.sqrt(level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_es**2)
                error_frac = error/o3
                o3_good = o3.where(mr>0.8).data
                ax.plot(o3*1e6, o3_z/1e3, '--', linewidth=0.1, color=color_spectro[spectro])
                #ax.plot(error_frac.where(mr>0.8).data*100, o3_z/1e3, '-.', linewidth=0.5, color=color_spectro[spectro])
                #ax.plot(error*1e6, o3_z/1e3, '-.', linewidth=0.5, color=color_spectro[spectro])
                ax.plot(o3_good*1e6, o3_z/1e3, linewidth=1, label=spectro,color=color_spectro[spectro])
                ax.set_title('Tb = ' + str(identifier_plot[i]), fontsize='small')
                ax.set_xlim(-0.5,12)
                ax.yaxis.set_major_locator(MultipleLocator(20))
                ax.yaxis.set_minor_locator(MultipleLocator(10))
                ax.xaxis.set_major_locator(MultipleLocator(5))
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                ax.grid(which='both',  axis='x', linewidth=0.5)
                ax.grid(which='minor', axis='y', linewidth=0.5)
                
                
                ax2.plot(mr, o3_z/1e3, linewidth=1, color=color_spectro[spectro], label=spectro)
                ax2.set_title('Tb = ' + str(identifier_plot[i]), fontsize='small')
                ax2.set_xlim(-0.1,1.2)
                ax2.yaxis.set_major_locator(MultipleLocator(20))
                ax2.yaxis.set_minor_locator(MultipleLocator(10))
                ax2.xaxis.set_major_locator(MultipleLocator(1))
                ax2.xaxis.set_minor_locator(MultipleLocator(0.1))
                ax2.grid(which='both',  axis='x', linewidth=0.5)
                ax2.grid(which='minor', axis='y', linewidth=0.5)

                #ax2.set_xticklabels([])
                if i==0 or i==5 or i==10:
                    ax.set_ylabel('alt [km]', fontsize='small')
                    ax2.set_ylabel('alt [km]', fontsize='small')
                else: 
                    ax.set_yticklabels([])
                    ax2.set_yticklabels([])

                if i in [10, 11, 12, 13, 14]:
                    ax.set_xlabel('$O_3$ [ppm]', fontsize='small')
                    ax2.set_xlabel('MR', fontsize='small') 
                else: 
                    ax.set_xticklabels([])
                    ax2.set_xticklabels([])

        ax.legend()
        fig.suptitle('$O_3$ VMR')
        figures2.append(fig)

        fig_mr.suptitle('Measurement response')
        figures2.append(fig_mr)

        #plt.plot(o3_apriori*1e6, o3_z/1e3)
        
        #fig, axs = plt.subplots(nrows=3, ncols=5, sharex=True, sharey=True)
        fig2 = plt.figure()
        fig2.subplots_adjust(hspace=0.4, wspace=0.4)
        for i in range(len(level2_data[spectro].observation)):
            for spectro in ['AC240','USRP-A']:
                ax = fig2.add_subplot(3,5,i+1)
                o3_diff = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_x - level2_data['U5303'].isel(observation=i, o3_lat=0, o3_lon=0).o3_x
                
                mr = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_mr
                mr_basis = level2_data['U5303'].isel(observation=i, o3_lat=0, o3_lon=0).o3_mr
                #o3_rel_diff = o3_diff / level2_data['U5303'].isel(observation=i, o3_lat=0, o3_lon=0).o3_x
                o3_z = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_z
                o3_p = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_p  
                o3_good_diff = o3_diff.where((mr.data >0.8) & (mr_basis.data>0.8))
                #o3_good_rel_diff = o3_rel_diff.where((mr.data >0.8) & (mr_basis.data>0.8))
                ax.plot(o3_diff*1e6, o3_z/1e3, '--', linewidth=0.2, color=color_spectro[spectro])
                ax.plot(o3_good_diff*1e6, o3_z/1e3, '-', linewidth=1, color=color_spectro[spectro], label=spectro)
                #ax.plot(o3_good_rel_diff*100, o3_z/1e3, '-', linewidth=1, color=color_spectro[spectro], label=spectro)
                ax.set_title('Tb = ' + str(identifier_plot[i]), fontsize='small')
                ax.set_xlim(-2,2)
                ax.set_ylim(10,70)
                ax.yaxis.set_major_locator(MultipleLocator(20))
                ax.yaxis.set_minor_locator(MultipleLocator(10))
                ax.xaxis.set_major_locator(MultipleLocator(2))
                ax.xaxis.set_minor_locator(MultipleLocator(0.2))
                ax.grid(which='minor',  axis='x', linewidth=0.2)
                ax.grid(which='major',  axis='x', linewidth=1)
                ax.grid(which='both', axis='y', linewidth=0.5)
                
                #ax2.set_xticklabels([])
                if i==0 or i==5 or i==10:
                    ax.set_ylabel('altitude [km]',fontsize='small')
                else: 
                    ax.set_yticklabels([])

                if i in [10, 11, 12, 13, 14]:
                    ax.set_xlabel('$O_3$ [ppm]',fontsize='small')
                else: 
                    ax.set_xticklabels([])
                    

        ax.legend()
        fig2.suptitle('$O_3$ VMR difference')
        figures2.append(fig2)

        fig4 = plt.figure()
        fig4.subplots_adjust(hspace=0.4, wspace=0.4)
        for spectro in spectro_lvl2:
            for i in range(len(level2_data[spectro].observation)):
                ax = fig4.add_subplot(3,5,i+1)
                obs_error = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_eo
                smoothing_error = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_es
                #error = np.sqrt(level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_es**2)
                o3_z = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_z
                o3_p = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_p  
                mr = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_mr

                ax.plot(obs_error*1e6, o3_z/1e3, '--', linewidth=0.5, color=color_spectro[spectro])
                ax.plot(obs_error.where(mr>0.8)*1e6, o3_z/1e3, '-.', linewidth=1, color=color_spectro[spectro], label='$e_o$, '+spectro)
                ax.plot(smoothing_error*1e6, o3_z/1e3, '--', linewidth=0.5, color=color_spectro[spectro])
                ax.plot(smoothing_error.where(mr>0.8)*1e6, o3_z/1e3, '-', linewidth=1, color=color_spectro[spectro], label='$e_s$,'+spectro)
                #ax.plot(error/1e3, o3_z/1e3, '-', linewidth=1, color=color_spectro[spectro], label=spectro)
                ax.set_title('Tb = ' + str(identifier_plot[i]), fontsize='small')
                #ax.set_xlim(-0.5,12)
                ax.yaxis.set_major_locator(MultipleLocator(20))
                ax.yaxis.set_minor_locator(MultipleLocator(10))
                ax.xaxis.set_major_locator(MultipleLocator(1))
                ax.xaxis.set_minor_locator(MultipleLocator(0.2))
                ax.grid(which='both',  axis='x', linewidth=0.5)
                ax.grid(which='minor', axis='y', linewidth=0.5)
                
                #ax2.set_xticklabels([])
                if i==0 or i==5 or i==10:
                    ax.set_ylabel('alt [km]',fontsize='small')
                else: 
                    ax.set_yticklabels([])

                if i in [10, 11, 12, 13, 14]:
                    ax.set_xlabel('$O_3$ [ppm]',fontsize='small')
                else: 
                    ax.set_xticklabels([])

        ax.legend()
        fig4.suptitle('Error')
        figures2.append(fig4)

        fig3 = plt.figure()
        fig3.subplots_adjust(hspace=0.4, wspace=0.4)
        for spectro in spectro_lvl2:
            for i in range(len(level2_data[spectro].observation)):
                ax = fig3.add_subplot(3,5,i+1)
                o3_offset = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_offset
                o3_fwhm = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_fwhm
                o3_z = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_z
                o3_p = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_p  
                #ax.plot(o3_offset/1e3, o3_z/1e3, '--', linewidth=0.5, color=color_spectro[spectro])
                ax.plot(o3_fwhm/1e3, o3_z/1e3, '-', linewidth=0.5, color=color_spectro[spectro], label=spectro)
                ax.set_title('Tb = ' + str(identifier_plot[i]), fontsize='small')
                #ax.set_xlim(-0.5,12)
                ax.yaxis.set_major_locator(MultipleLocator(20))
                ax.yaxis.set_minor_locator(MultipleLocator(10))
                ax.xaxis.set_major_locator(MultipleLocator(10))
                ax.xaxis.set_minor_locator(MultipleLocator(5))
                ax.grid(which='both',  axis='x', linewidth=0.5)
                ax.grid(which='minor', axis='y', linewidth=0.5)
                
                #ax2.set_xticklabels([])
                if i==0 or i==5 or i==10:
                    ax.set_ylabel('alt [km]',fontsize='small')
                else: 
                    ax.set_yticklabels([])

                if i in [10, 11, 12, 13, 14]:
                    ax.set_xlabel('FWHM [km]',fontsize='small')
                else: 
                    ax.set_xticklabels([])

        ax.legend()
        fig3.suptitle('FWHM')
        figures2.append(fig3)
        save_single_pdf(basename_lvl2+'bias_o3_feb'+'_'+integration_strategy+'.pdf',figures2)
# %%    