#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:58:02 2020

@author: eric

Module with function dealing with mopi5 data at every level

"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

from GROSOM_library import tropospheric_correction
from utils_GROSOM import save_single_pdf

import matplotlib.dates as mdates
import datetime

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.lines import Line2D

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]})

color_spectro = {'AC240':'tab:orange', 'USRP-A':'tab:green', 'U5303':'tab:blue', 'AC240_unbiased':'tab:red'}
F0 = 110.836e9

def return_bad_channels_mopi5(number_of_channel, date, spectro):
    '''
    common function for the 2 classes... maybe a better way to do it.
    '''
    if spectro == 'USRP-A':
        bad_channels = np.hstack((np.arange(0,1023), 8092, np.arange(number_of_channel-1023,number_of_channel)))
    elif spectro == 'U5303':
        #bad_channels = np.where(intermediate_frequency > 1000)
        bad_channels = np.hstack((np.arange(0,63), np.arange(10500,number_of_channel)))
    elif spectro == 'AC240':
        bad_channels = np.hstack((np.arange(0,63), np.arange(number_of_channel-63,number_of_channel)))
    else:
        ValueError('Spectrometer unknown !')
    
    return bad_channels

def compare_Tb_mopi5(self, ds_dict, calibration_cycle):
    fig, axs = plt.subplots(1,1,sharex=True)
    for s in self.spectrometers:
        mask = ds_dict[s].good_channels[calibration_cycle].data
        mask[mask==0]=np.nan
        axs.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,ds_dict[s].Tb[calibration_cycle]*mask, lw=0.5, label=s)
        axs.set_xlim(110.25, 111.4)
        #axs.set_ylim(0,250)
        #ax].set_ylim(np.median(ds_dict[s].Tb[calibration_cycle].data)-10,np.median(ds_dict[s].Tb[calibration_cycle].data)+15)
        axs.set_xlabel("f [GHz]")
        axs.set_ylabel(r"$T_B$ [K]")
        axs.set_title("Tb")
        axs.grid()
        axs.legend(fontsize='xx-small')
        #ax3.legend()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    return fig

def correct_troposphere(calibration, spectrometers, dim, method='Ingold_v1', use_basis='AC240', skip_ch = [1000,1000], num_of_ch=500):
    '''
    generic tropospheric correction for GROSOM
    '''

    # little trick to work with other dimensions than time
    tb_corr_da = dict()
    mean_opacity = dict()
    mean_transmittance = dict()

    for s in spectrometers:
        if dim != 'time':
            calibration.integrated_data[s] = calibration.integrated_data[s].swap_dims({dim:'time'})
            calibration.integrated_meteo[s] = calibration.integrated_meteo[s].swap_dims({dim:'time'})

        tb_corr_da[s] = xr.DataArray(None,
            coords = [calibration.integrated_data[s].coords['time'], calibration.integrated_data[s].coords['channel_idx']])
        mean_opacity[s] = xr.DataArray(None,
            coords = [calibration.integrated_data[s].coords['time']])
        mean_transmittance[s] = xr.DataArray(None,
            coords = [calibration.integrated_data[s].coords['time']])
    
    # first run with the spectrometer to use as basis (just initiating the corr_func)
    for i, t in enumerate(calibration.integrated_data[s].coords['time']):
        s = use_basis
        try:
            temperature = calibration.integrated_meteo[s].air_temperature.isel(time=i).data.item()
        except:
            print('No ground temperature found, taking default')
            temperature = 293
        if method == 'Ingold_v1':
            delta_T = 10.4
            tropospheric_temperature = temperature - delta_T
        elif method == 'Ingold_v2':
            raise NotImplementedError
        
        bad_channels = calibration.integrated_data[s].good_channels.isel(time=i).data
        bad_channels[bad_channels==0]=np.nan
        clean_tb = calibration.integrated_data[s].Tb.isel(time=i).data * bad_channels
        
        corr_func = tropospheric_correction(
            f = calibration.integrated_data[s].frequencies.isel(time=i).data,
            y = clean_tb,
            t_trop=tropospheric_temperature,
            use_wings='both',
            skip=skip_ch,
            num_el=num_of_ch)
    
        # Run correction on all spectrometers
        for s in spectrometers:        
        #for i, t in enumerate(calibration.integrated_data[s].coords['time']):
            tb_corr, opac, transm = corr_func(
                f = calibration.integrated_data[s].frequencies.isel(time=i).data,
                y = calibration.integrated_data[s].Tb.isel(time=i).data
                )
            tb_corr_da[s][i,:] =  tb_corr
            mean_transmittance[s][i] = transm
            mean_opacity[s][i] = opac

    for s in spectrometers:    
        calibration.integrated_data[s] = calibration.integrated_data[s].assign(Tb_corr = tb_corr_da[s])
        calibration.integrated_data[s] = calibration.integrated_data[s].assign(tropospheric_transmittance = mean_transmittance[s])
        calibration.integrated_data[s] = calibration.integrated_data[s].assign(tropospheric_opacity = mean_opacity[s])
        
        if dim != 'time':
            calibration.integrated_data[s] = calibration.integrated_data[s].swap_dims({'time':dim})
            calibration.integrated_meteo[s] = calibration.integrated_meteo[s].swap_dims({'time':dim})
    
    return calibration.integrated_data


def correct_troposphere_interpolated(calibration, spectrometers, dim, method='Ingold_v1', use_basis='AC240', skip_ch = [1000,1000], num_of_ch=500):
    '''
    generic tropospheric correction for GROSOM
    '''
    
    # little trick to work with other dimensions than time
    tb_corr_da = dict()
    mean_opacity = dict()
    mean_transmittance = dict()
    tb_corr_da_interp = dict()


    for s in spectrometers:
        tb_corr_da[s] = xr.DataArray(None,
            coords = [calibration.integrated_data[s].coords[dim], calibration.integrated_data[s].coords['channel_idx']])
        mean_opacity[s] = xr.DataArray(None,
            coords = [calibration.integrated_data[s].coords[dim]])
        mean_transmittance[s] = xr.DataArray(None,
            coords = [calibration.integrated_data[s].coords[dim]])
        
        tb_corr_da_interp[s] = xr.DataArray(None,
            coords = [calibration.integrated_data[s].coords[dim], calibration.integrated_data[s].coords['bin_freq_interp']])

        if dim != 'time':
            calibration.integrated_data[s] = calibration.integrated_data[s].swap_dims({dim:'time'})
            calibration.integrated_meteo[s] = calibration.integrated_meteo[s].swap_dims({dim:'time'})
    
    # first run with the spectrometer to use as basis (just initiating the corr_func)
    for i, t in enumerate(calibration.integrated_data[s].coords[dim]):
        s = use_basis
        try:
            temperature = calibration.integrated_meteo[s].air_temperature.isel(time=i).data.item()
        except:
            print('No ground temperature found, taking default')
            temperature = 293
        if method == 'Ingold_v1':
            delta_T = 10.4
            tropospheric_temperature = temperature - delta_T
        elif method == 'Ingold_v2':
            raise NotImplementedError
        
        bad_channels = calibration.integrated_data[s].good_channels.isel(time=i).data
        bad_channels[bad_channels==0]=np.nan
        clean_tb = calibration.integrated_data[s].Tb.isel(time=i).data * bad_channels
        
        clean_tb_interp = calibration.integrated_data[s].interpolated_Tb.isel(time=i).data
        
        corr_func = tropospheric_correction(
            f = calibration.integrated_data[s].frequencies.isel(time=i).data,
            y = clean_tb,
            t_trop=tropospheric_temperature,
            use_wings='both',
            skip=skip_ch,
            num_el=num_of_ch)
    
        # Run correction on all spectrometers
        for s in spectrometers:        
        #for i, t in enumerate(calibration.integrated_data[s].coords['time']):
            tb_corr_interp, opac, transm = corr_func(
                f = calibration.integrated_data[s].bin_freq_interp.data,
                y = clean_tb_interp
                )
            tb_corr_da_interp[s][i,:] =  tb_corr_interp
            # mean_transmittance[s][i] = transm
            # mean_opacity[s][i] = opac

    for s in spectrometers: 
        if dim != 'time':
            calibration.integrated_data[s] = calibration.integrated_data[s].swap_dims({'time':dim})
            #calibration.integrated_meteo[s] = calibration.integrated_data[s].swap_dims({'time':dim})
        #if dim != 'time':
            # tb_corr_da[s] = tb_corr_da[s].swap_dims({'time':dim})
            # mean_transmittance[s] = mean_transmittance[s].swap_dims({'time':dim})
            # mean_opacity[s] = mean_opacity[s].swap_dims({'time':dim})

        calibration.integrated_data[s] = calibration.integrated_data[s].assign(Tb_corr_interp = tb_corr_da_interp[s])
        #calibration.integrated_data[s] = calibration.integrated_data[s].assign(tropospheric_transmittance_interp = mean_transmittance[s])
        #calibration.integrated_data[s] = calibration.integrated_data[s].assign(tropospheric_opacity_interp = mean_opacity[s])
        
    return calibration.integrated_data

def compare_spectra_mopi5_new(cal_int_obj, ds_dict, spectrometers, calibration_cycle=0, title='', corr_band=False):
    #fig, axs = plt.subplots(2,2,sharex=True)
    fig = plt.figure()
    ax1 = fig.add_subplot(2,2, (1,2))
    ax2 = fig.add_subplot(2,2, (3,4))
    ax3 = ax1.inset_axes([0.1, 0.6, 0.15, 0.35])
    ax4 = ax2.inset_axes([0.1, 0.6, 0.15, 0.35])
    for s in spectrometers:
        mask = ds_dict[s].good_channels[calibration_cycle].data
        mask[mask==0]=np.nan
        ax1.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,
                 ds_dict[s].Tb[calibration_cycle].data*mask, lw=0.5, label=s)
        if corr_band:
            if s=='AC240':
                freq_left = [ds_dict[s].frequencies[calibration_cycle].data[2459]/1e9, ds_dict[s].frequencies[calibration_cycle].data[3259]/1e9]
                freq_right = [ds_dict[s].frequencies[calibration_cycle].data[-818]/1e9, ds_dict[s].frequencies[calibration_cycle].data[-1618]/1e9]
                c='orange'
            elif s == 'U5303':
                freq_left = [ds_dict[s].frequencies[calibration_cycle].data[1538]/1e9, ds_dict[s].frequencies[calibration_cycle].data[2038]/1e9]
                freq_right = [ds_dict[s].frequencies[calibration_cycle].data[-6654]/1e9, ds_dict[s].frequencies[calibration_cycle].data[-7154]/1e9]
                c='blue'
            ax1.axvspan(freq_left[0], freq_left[1], alpha=0.5, color=c)
            ax1.axvspan(freq_right[0], freq_right[1], alpha=0.5, color=c)

        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.set_title(title)
        ax1.grid()
        #ax1.legend(fontsize='xx-small')
        # axs[1][0].plot(ds_dict[s].frequencies.data/1e9,
        #          ds_dict[s].Tb_corr_old[id].data*mask,
        #          lw=0.5, label=s)
        # #axs[1][0].set_xlim(110.3, 111.4)
        # axs[1][0].set_ylim(0,30)
        # axs[1][0].set_xlabel("f [GHz]")
        # axs[1][0].set_ylabel(r"$T_B$ [K]")
        # axs[1][0].set_title("Tb_corr_old")
        # axs[1][0].grid()
        ax2.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,
                 ds_dict[s].Tb_corr[calibration_cycle].data*mask,
                 lw=0.5, label=s)
        ax2.set_xlim(110.25, 111.4)
        ax2.set_ylim(0,30)
        ax2.set_xlabel("f [GHz]")
        ax2.set_ylabel(r"$T_B$ [K]")
        ax2.set_title("Tb_corr")
        ax2.legend(fontsize='xx-small',loc=1)
        ax2.grid()
        ax3.plot(
            (ds_dict[s].frequencies[calibration_cycle].data-cal_int_obj.observation_frequency)/1e6,
            ds_dict[s].Tb[calibration_cycle].data*mask, lw=0.2
        )
        ax3.set_xlim(-10, 10)
        #ax3.set_xlabel('[MHz]')
        ymax = np.nanmax(ds_dict[s].Tb[calibration_cycle].data*mask)
        if not np.isnan(ymax):
            ax3.set_ylim(ymax-5,ymax+0.5)
        ax3.grid()
        ax4.plot(
            (ds_dict[s].frequencies[calibration_cycle].data-cal_int_obj.observation_frequency)/1e6,
            ds_dict[s].Tb_corr[calibration_cycle].data*mask, lw=0.2
        )
        ax4.set_xlim(-10, 10)
        #ax4.set_xlabel('[MHz]')
        ymax = np.nanmax(ds_dict[s].Tb_corr[calibration_cycle].data*mask)
        if not np.isnan(ymax):
            ax4.set_ylim(ymax-5,ymax+0.5)
        ax4.grid()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def plot_time_min_comp(cal_int_obj, title='Mean of integration time'):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    symbol = ['x','o','.']
    for i, s in enumerate(cal_int_obj.spectrometers):
        ax1.plot(cal_int_obj.integrated_data[s].time_min, symbol[i], label=s)
    ax1.legend()

def compare_spectra_only_mopi5(cal_int_obj, ds_dict, spectrometers, calibration_cycle=0, title='', corr_band=[]):
    #print(corr_band)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    col = 'k'
   # ax2 = ax1.inset_axes([0.1, 0.5, 0.3, 0.45])
    for s in spectrometers:
        mask = ds_dict[s].good_channels[calibration_cycle].data
        mask[mask==0]=np.nan
        ax1.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,
                 ds_dict[s].Tb[calibration_cycle].data*mask, color=col, lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.3)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[calibration_cycle].data)-10,np.median(ds_dict[s].Tb[id].data)+calibration_cycle)
        ax1.set_xlabel("Frequency [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(2))
        ax1.set_title(title)
        
       #s ax1.legend(fontsize='small',loc=1)
        if corr_band:
            freq_left=ds_dict[s].frequencies[calibration_cycle].where(
                abs(ds_dict[s].frequencies[calibration_cycle]-corr_band[s][0])<corr_band[s][1],drop=True)
            freq_right=ds_dict[s].frequencies[calibration_cycle].where(
                abs(ds_dict[s].frequencies[calibration_cycle]-corr_band[s][2])<corr_band[s][3],drop=True)
            if (s=='AC240') or (s=='U5303'):
                b = ds_dict[s].continuum_value_line_center[calibration_cycle].data - ds_dict[s].slope_indiv[calibration_cycle].data*cal_int_obj.observation_frequency
                ax1.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,
                ds_dict[s].slope_indiv[calibration_cycle].data*ds_dict[s].frequencies[calibration_cycle].data+b, color=col, lw=0.5)
                ax1.axvspan(freq_left[0]/1e9, freq_left[-1]/1e9 , alpha=0.25, color=col)
                ax1.axvspan(freq_right[0]/1e9, freq_right[-1]/1e9, alpha=0.25, color=col)
        

        #ax3.legend()
        #ax2.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,ds_dict[s].Tb[calibration_cycle].data*mask, lw=0.2)
        # ax2.plot(
        #     (ds_dict[s].frequencies[calibration_cycle].data-cal_int_obj.observation_frequency)/1e6,
        #     ds_dict[s].Tb[calibration_cycle].data*mask, lw=0.2
        # )
        # ax2.yaxis.set_major_locator(MultipleLocator(2))
        # ax2.yaxis.set_minor_locator(MultipleLocator(1))
        # ax2.set_xlim(-10, 10)
        # ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        # ymax = np.nanmax(ds_dict[s].Tb[calibration_cycle].data*mask)
        # if not np.isnan(ymax):
        #     ax2.set_ylim(ymax-3,ymax+1)
        
        #ax2.legend(fontsize='xx-small')
    ax1.grid()
    # ax2.grid(which='both')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def compare_binned_spectra_only_mopi5(cal_int_obj, ds_dict, calibration_cycle=0, use_basis='U5303', title=''):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = ax1.inset_axes([0.1, 0.6, 0.15, 0.35])
    ax3 = fig.add_subplot(212)
    for s in cal_int_obj.spectrometers:
        ax1.plot(ds_dict[s].bin_freq.data/1e9, ds_dict[s].binned_Tb[calibration_cycle], lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.set_title(title)
        ax1.grid()
        #ax1.legend(fontsize='xx-small')
        #ax3.legend()
        #ax2.plot(ds_dict[s].bin_freq.data/1e9, ds_dict[s].binned_Tb[calibration_cycle], lw=0.2)
        ax2.plot(
            (ds_dict[s].bin_freq.data-cal_int_obj.observation_frequency)/1e6,
            ds_dict[s].binned_Tb[calibration_cycle], lw=0.2
        )
        ax2.set_xlim(-10, 10)
        #ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        ymax = np.nanmax(ds_dict[s].binned_Tb[calibration_cycle])
        if not np.isnan(ymax):
            ax2.set_ylim(ymax-4,ymax+0.5)
        ax2.grid()

        Tb_diff = ds_dict[use_basis].binned_Tb[calibration_cycle].data-ds_dict[s].binned_Tb[calibration_cycle].data
        
        ax3.plot(ds_dict[s].bin_freq.data/1e9, Tb_diff, lw=0.5, label=s)
        ax3.set_title('Tb difference with: '+use_basis)
        ax3.set_ylim(-0.5,1.5)
        ax3.legend(fontsize='xx-small')
        ax3.grid()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def plot_ts_mopi5(calibration, title):
    import matplotlib.dates as mdates
    fig = plt.figure()
    ax = fig.subplots(nrows=3, ncols=1, sharex=True)
    #ax2 = fig.add_subplot(412, sharex=True)
    ax31 = ax[1]
    ax32 = ax31.twinx()
    ax41 = ax[2]
    ax42 = ax41.twinx()
    #ax2 = fig.add_subplot(3,1,3)
    meteo_lw = 0.5
    marker = ['x','*','.']
    for a, spectro in enumerate(calibration.spectrometers):
        #dates = [pd.to_datetime(d) for d in calibration.calibrated_data[spectro].time.data]
        dates =calibration.calibrated_data[spectro].time.data
        ax[0].scatter(dates, calibration.calibrated_data[spectro].mean_Tb, marker=marker[a], s=2, label=spectro)
        ax[0].set_xlim(dates[0],dates[-1])
        #ax[1].plot(calibration.calibrated_data[s].time, calibration.calibrated_data[s].TSys,'.',markersize=4,label=s)
        #ax2.plot(calibration.calibration_flags[s].time, calibration.calibration_flags[s].sum(dim='flags').to_array()[0,:],'.',markersize=4,label=s)
        #ax.legend()
    ax[0].set_ylim(70,260)
    ax[0].legend(fontsize='xx-small')
    s = 'AC240'

    ax31.plot(calibration.meteo_data[s].time, calibration.meteo_data[s].air_temperature-273.15,'r-', lw=meteo_lw) 
    ax32.plot(calibration.meteo_data[s].time, calibration.meteo_data[s].air_pressure,'k-', lw=meteo_lw)  
    ax31.set_ylabel('$T_{air}$ [$\degree C$]', color='r')
    ax31.set_ylim(-6, 22)
    ax31.tick_params(axis='y', labelcolor='r')

    ax32.set_ylabel('$P_{air}$ [hPa]', color='k')
    ax32.set_ylim(920, 970)

    ax42.plot(calibration.meteo_data[s].time, calibration.meteo_data[s].relative_humidity,'b-', lw=meteo_lw)
    ax41.bar(calibration.meteo_data[s].time.data, calibration.meteo_data[s].precipitation.data, width=0.01, color='k')
    

    ax42.set_ylabel('$RH$ [%]', color='b')
    ax42.set_ylim(20, 100)
    ax42.yaxis.set_major_locator(MultipleLocator(20))
    ax42.tick_params(axis='y', labelcolor='b')
    

    ax41.set_ylabel('Prec [mm]', color='k')
    ax41.tick_params(axis='y', labelcolor='k')
    ax41.set_ylim(0, 0.2)
    ax41.yaxis.set_major_locator(MultipleLocator(0.05))
    #ax41.yaxis.set_minor_locator(MultipleLocator(10))
    
    ax[0].set_ylabel('$T_B$ [K]')

    days = mdates.DayLocator() 
    hours = mdates.HourLocator() 
    
    date_form = mdates.DateFormatter("%m-%d")

    # ax.yaxis.set_major_locator(MultipleLocator(20))
    # ax.yaxis.set_minor_locator(MultipleLocator(10))
    # ax.xaxis.set_major_locator(MultipleLocator(10))
    # ax.xaxis.set_minor_locator(MultipleLocator(5))
    # ax.grid(which='both',  axis='x', linewidth=0.5)
    # ax.grid(which='minor', axis='y', linewidth=0.5)
    # #ax[1].xaxis.set_major_formatter(date_form)

    # Ensure a major tick for each week using (interval=1) 
    #ax[1].xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
    
    ax[0].yaxis.set_major_locator(MultipleLocator(50))
    ax[0].yaxis.set_minor_locator(MultipleLocator(10))
    ax[0].xaxis.set_major_locator(days)
    ax[0].xaxis.set_minor_locator(hours)   

    ax[1].xaxis.set_major_locator(days)
    ax[1].xaxis.set_minor_locator(hours)
    #plt.suptitle('Mean T_B and T_sys')
    for i in range(len(ax)):
        ax[i].grid(which='major')
        ax[i].grid(which='minor', linewidth=0.2)
    fig.tight_layout(rect=[0, 0.01, 1, 1])
    #fig.savefig('/scratch/MOPI5/Level1/time_series_2019_04.png', dpi=600, facecolor='w', edgecolor='w')
    
    return fig


def plot_ts_mopi5_Feb_paper(calibration, title):
    fig = plt.figure()
    ax = fig.subplots(nrows=3, ncols=1, sharex=True)
    #ax2 = fig.add_subplot(412, sharex=True)
    ax31 = ax[1]
    ax32 = ax31.twinx()
    ax41 = ax[2]
    ax42 = ax41.twinx()
    #ax2 = fig.add_subplot(3,1,3)
    meteo_lw = 0.75
    marker = ['x','*','.']

    start = datetime.date(2019,1,30)
    stop = datetime.date(2019,2,23)

    for a, spectro in enumerate(calibration.spectrometers):
        #dates = [pd.to_datetime(d) for d in calibration.calibrated_data[spectro].time.data]
        dates =calibration.calibrated_data[spectro].time.data
        xr.plot.scatter(calibration.calibrated_data[spectro], x= 'time', y='mean_Tb', ax=ax[0], marker=marker[a], s=2, label=spectro)
        # ax[0].scatter(dates, calibration.calibrated_data[spectro].mean_Tb, marker=marker[a], s=2, label=spectro)
        #ax[1].plot(calibration.calibrated_data[s].time, calibration.calibrated_data[s].TSys,'.',markersize=4,label=s)
        #ax2.plot(calibration.calibration_flags[s].time, calibration.calibration_flags[s].sum(dim='flags').to_array()[0,:],'.',markersize=4,label=s)
        #ax.legend()
    ax[0].set_xlim([start,stop])
    ax[0].set_ylim(60,260)
    ax[0].set_ylabel('$T_B$ [K]')
    ax[0].set_xlabel('')

    ax[0].legend(fontsize='small', loc='upper left', bbox_to_anchor=(0.75, 0.99))
    s = 'AC240'

    #ts1 = calibration.meteo_data[s].where(calibration.meteo_data[s].time < pd.to_datetime('2019-02-02'), drop=True)
    ts1 = calibration.meteo_data[s].sel(time=slice('2019-01-30','2019-02-22'), drop=True)

    # ts2 = calibration.meteo_data[s].sel(time=slice('2019-02-03','2019-02-10'),drop=True)

    # #ts2 = calibration.meteo_data[s].where(calibration.meteo_data[s].time > pd.to_datetime('2019-02-02'), drop=True)
    # #ts2 = ts2.where(calibration.meteo_data[s].time < pd.to_datetime('2019-02-11'), drop=True)
    
    # #ts3 = calibration.meteo_data[s].where(calibration.meteo_data[s].time > pd.to_datetime('2019-02-11'), drop=True)
    # ts3 = calibration.meteo_data[s].sel(time=slice('2019-02-11','2019-02-22'),drop=True)
    
    (ts1.air_temperature-273.15).resample(time='1H', skipna=True).mean().plot.line(
        'red',
        ax=ax31,
        xticks=[],
        lw=meteo_lw
    )
    ts1.air_pressure.resample(time='1H', skipna=True).mean().plot.line(
        'k',
        ax=ax32,
        xticks=[],
        lw=meteo_lw
    )
    

    # ax31.plot(ts1.time, ts1.air_temperature-273.15,'r-', lw=meteo_lw) 
    # ax31.plot(ts2.time, ts2.air_temperature-273.15,'r-', lw=meteo_lw) 
    # ax31.plot(ts3.time, ts3.air_temperature-273.15,'r-', lw=meteo_lw) 

    # ax32.plot(ts1.time, ts1.air_pressure,'k-', lw=meteo_lw)  
    # ax32.plot(ts2.time, ts2.air_pressure,'k-', lw=meteo_lw)  
    # ax32.plot(ts3.time, ts3.air_pressure,'k-', lw=meteo_lw)  
    ax31.set_ylabel('$T_{air}$ [$\degree C$]', color='r')
    ax31.set_ylim(-10, 20)
    ax31.tick_params(axis='y', labelcolor='r')

    ax32.set_ylabel('$P_{air}$ [hPa]', color='k')
    ax32.set_ylim(910, 970)

    ax31.set_xlabel('')
    ax32.set_xlabel('')

    ax31.yaxis.set_major_locator(MultipleLocator(5))
    ax32.yaxis.set_major_locator(MultipleLocator(20))
    # ax42.plot(ts1.time, ts1.relative_humidity,'b-', lw=meteo_lw)
    # ax42.plot(ts2.time, ts2.relative_humidity,'b-', lw=meteo_lw)
    # ax42.plot(ts3.time, ts3.relative_humidity,'b-', lw=meteo_lw)
    
    #rain_serie = calibration.meteo_data[s].precipitation.resample(time='2H', skipna=True).sum().to_series()
    #rain_serie.plot.bar(ax=ax41, x='time', width=0.1, color='k',xticks=[])

    ts1.relative_humidity.resample(time='1H', skipna=True).mean().plot.line(
        'b-',
        ax=ax42,
        xticks=[],
        lw=meteo_lw
    )

    rain = calibration.meteo_data[s].precipitation.resample(time='1H', skipna=True).sum()
    ax41.bar(rain.time.data, rain.data, width=0.08, color='k')
    
    ax42.set_ylabel('$RH$ [%]', color='b')
    ax42.tick_params(axis='y', labelcolor='b')
    ax42.set_ylim(20, 100)
    ax42.yaxis.set_major_locator(MultipleLocator(20))

    ax41.set_ylabel('Prec [mm]', color='k')
    ax41.tick_params(axis='y', labelcolor='k')
    ax41.set_ylim(0, 2)
    ax41.yaxis.set_major_locator(MultipleLocator(0.5))
    #ax41.yaxis.set_minor_locator(MultipleLocator(10))
    

    midmonth = mdates.DayLocator(7) 
    days = mdates.DayLocator() 
    hours = mdates.HourLocator() 
    
    #ax[1].xaxis.set_major_locator(days)
    date_form = mdates.DateFormatter("%m-%d")

    # ax.yaxis.set_major_locator(MultipleLocator(20))
    # ax.yaxis.set_minor_locator(MultipleLocator(10))
    # ax.xaxis.set_major_locator(MultipleLocator(10))
    # ax.xaxis.set_minor_locator(MultipleLocator(5))
    # ax.grid(which='both',  axis='x', linewidth=0.5)
    # ax.grid(which='minor', axis='y', linewidth=0.5)
    # #ax[1].xaxis.set_major_formatter(date_form)

    # Ensure a major tick for each week using (interval=1) 
    ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=5))
    
    ax[0].yaxis.set_major_locator(MultipleLocator(50))
    ax[0].yaxis.set_minor_locator(MultipleLocator(10))
    
    #ax[1].xaxis.set_major_locator(midmonth)
    ax[1].xaxis.set_minor_locator(days)
    #plt.suptitle('Mean T_B and T_sys')
    for i in range(len(ax)):
        ax[i].grid(which='major')
        ax[i].grid(which='minor', linewidth=0.2)
    fig.tight_layout(rect=[0, 0.01, 1, 1])
    #fig.savefig('/scratch/MOPI5/Level1/time_series_2019_02_22.png', dpi=600, facecolor='w', edgecolor='w')

    return fig


def plot_ts_mopi5_Feb(calibration, title):
    fig = plt.figure()
    ax = fig.subplots(nrows=3, ncols=1, sharex=True)
    #ax2 = fig.add_subplot(412, sharex=True)
    ax31 = ax[1]
    ax32 = ax31.twinx()
    ax41 = ax[2]
    ax42 = ax41.twinx()
    #ax2 = fig.add_subplot(3,1,3)
    meteo_lw = 0.5
    marker = ['x','*','.']

    start = datetime.date(2019,1,30)
    stop = datetime.date(2019,2,23)

    for a, spectro in enumerate(calibration.spectrometers):
        #dates = [pd.to_datetime(d) for d in calibration.calibrated_data[spectro].time.data]
        dates =calibration.calibrated_data[spectro].time.data
        ax[0].scatter(dates, calibration.calibrated_data[spectro].mean_Tb, marker=marker[a], s=2, label=spectro)
        #ax[1].plot(calibration.calibrated_data[s].time, calibration.calibrated_data[s].TSys,'.',markersize=4,label=s)
        #ax2.plot(calibration.calibration_flags[s].time, calibration.calibration_flags[s].sum(dim='flags').to_array()[0,:],'.',markersize=4,label=s)
        #ax.legend()
    ax[0].set_xlim([start,stop])
    ax[0].set_ylim(60,260)

    ax[0].legend(fontsize='small', loc='upper left', bbox_to_anchor=(0.75, 0.99))
    s = 'AC240'

    #ts1 = calibration.meteo_data[s].where(calibration.meteo_data[s].time < pd.to_datetime('2019-02-02'), drop=True)
    ts1 = calibration.meteo_data[s].sel(time=slice('2019-01-30','2019-02-02'),drop=True)

    ts2 = calibration.meteo_data[s].sel(time=slice('2019-02-03','2019-02-10'),drop=True)

    #ts2 = calibration.meteo_data[s].where(calibration.meteo_data[s].time > pd.to_datetime('2019-02-02'), drop=True)
    #ts2 = ts2.where(calibration.meteo_data[s].time < pd.to_datetime('2019-02-11'), drop=True)
    
    #ts3 = calibration.meteo_data[s].where(calibration.meteo_data[s].time > pd.to_datetime('2019-02-11'), drop=True)
    ts3 = calibration.meteo_data[s].sel(time=slice('2019-02-11','2019-02-22'),drop=True)
    
    ax31.plot(ts1.time, ts1.air_temperature-273.15,'r-', lw=meteo_lw) 
    ax31.plot(ts2.time, ts2.air_temperature-273.15,'r-', lw=meteo_lw) 
    ax31.plot(ts3.time, ts3.air_temperature-273.15,'r-', lw=meteo_lw) 

    ax32.plot(ts1.time, ts1.air_pressure,'k-', lw=meteo_lw)  
    ax32.plot(ts2.time, ts2.air_pressure,'k-', lw=meteo_lw)  
    ax32.plot(ts3.time, ts3.air_pressure,'k-', lw=meteo_lw)  
    ax31.set_ylabel('$T_{air}$ [$\degree C$]', color='r')
    ax31.set_ylim(-6, 22)
    ax31.tick_params(axis='y', labelcolor='r')

    ax32.set_ylabel('$P_{air}$ [hPa]', color='k')
    ax32.set_ylim(920, 970)

    ax42.plot(ts1.time, ts1.relative_humidity,'b-', lw=meteo_lw)
    ax42.plot(ts2.time, ts2.relative_humidity,'b-', lw=meteo_lw)
    ax42.plot(ts3.time, ts3.relative_humidity,'b-', lw=meteo_lw)
    
    ax41.bar(calibration.meteo_data[s].time.data, calibration.meteo_data[s].precipitation.data, width=0.01, color='k')
    
    ax42.set_ylabel('$RH$ [%]', color='b')
    ax42.tick_params(axis='y', labelcolor='b')
    ax42.set_ylim(20, 100)
    ax42.yaxis.set_major_locator(MultipleLocator(20))

    ax41.set_ylabel('Prec [mm]', color='k')
    ax41.tick_params(axis='y', labelcolor='k')
    ax41.set_ylim(0, 0.2)
    ax41.yaxis.set_major_locator(MultipleLocator(0.05))
    #ax41.yaxis.set_minor_locator(MultipleLocator(10))
    
    ax[0].set_ylabel('$T_B$ [K]')

    midmonth = mdates.DayLocator(7) 
    days = mdates.DayLocator() 
    hours = mdates.HourLocator() 
    
    #ax[1].xaxis.set_major_locator(days)
    date_form = mdates.DateFormatter("%m-%d")

    # ax.yaxis.set_major_locator(MultipleLocator(20))
    # ax.yaxis.set_minor_locator(MultipleLocator(10))
    # ax.xaxis.set_major_locator(MultipleLocator(10))
    # ax.xaxis.set_minor_locator(MultipleLocator(5))
    # ax.grid(which='both',  axis='x', linewidth=0.5)
    # ax.grid(which='minor', axis='y', linewidth=0.5)
    # #ax[1].xaxis.set_major_formatter(date_form)

    # Ensure a major tick for each week using (interval=1) 
    ax[1].xaxis.set_major_locator(mdates.DayLocator(interval=5))
    
    ax[0].yaxis.set_major_locator(MultipleLocator(50))
    ax[0].yaxis.set_minor_locator(MultipleLocator(10))
    
    #ax[1].xaxis.set_major_locator(midmonth)
    ax[1].xaxis.set_minor_locator(days)
    #plt.suptitle('Mean T_B and T_sys')
    for i in range(len(ax)):
        ax[i].grid(which='major')
        ax[i].grid(which='minor', linewidth=0.2)
    fig.tight_layout(rect=[0, 0.01, 1, 1])
    #fig.savefig('/scratch/MOPI5/Level1/time_series_2019_02_22.png', dpi=600, facecolor='w', edgecolor='w')

    return fig


def compare_spectra_diff_mopi5(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', title=''):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = ax1.inset_axes([0.1, 0.6, 0.15, 0.35])
    ax3 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].Tb[calibration_cycle].where(ds_dict[use_basis].good_channels[calibration_cycle]==1).data
    clean_f = ds_dict[use_basis].frequencies[calibration_cycle].data
    for s in spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb[calibration_cycle].data
        Tb_diff = clean_Tb-Tb
        ax1.plot(clean_f/1e9, ds_dict[s].interpolated_Tb[calibration_cycle].data, lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.set_title(title)
        ax1.grid()
        #ax1.legend(fontsize='xx-small')
        #ax3.legend()
        ax2.plot((clean_f-cal_int_obj.observation_frequency)/1e6, ds_dict[s].interpolated_Tb[calibration_cycle].data, lw=0.2)
        ax2.set_xlim(-10, 10)
        #ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        ymax = np.nanmax(ds_dict[s].interpolated_Tb[calibration_cycle].data)
        if not np.isnan(ymax):
            ax2.set_ylim(ymax-4,ymax+0.5)
        ax2.grid()
        ax3.plot(clean_f/1e9, Tb_diff, lw=0.5, label=s)
        ax3.set_title('Interpolated Tb difference with: '+use_basis)
        ax3.set_ylim(-0.5,1.5)
        ax3.grid()
        ax3.legend(fontsize='xx-small')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def compare_spectra_binned_interp_mopi5(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', title=''):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = ax1.inset_axes([0.1, 0.6, 0.15, 0.35])
    ax3 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].binned_Tb[calibration_cycle].data
    clean_f = ds_dict[use_basis].bin_freq.data
    for s in spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb[calibration_cycle].data
        Tb_diff = Tb-clean_Tb
        ax1.plot(clean_f/1e9, ds_dict[s].interpolated_Tb[calibration_cycle].data, lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.set_title(title)
        ax1.grid()
        #ax1.legend(fontsize='xx-small')
        #ax3.legend()
        ax2.plot((clean_f-cal_int_obj.observation_frequency)/1e6, ds_dict[s].interpolated_Tb[calibration_cycle].data, lw=0.2)
        ax2.set_xlim(-10, 10)
        #ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        ymax = np.nanmax(ds_dict[s].interpolated_Tb[calibration_cycle].data)
        if not np.isnan(ymax):
            ax2.set_ylim(ymax-4,ymax+0.5)
        ax2.grid()
        ax3.plot(clean_f/1e9, Tb_diff, lw=0.5, label=s)
        ax3.set_title('Binned and interpolated Tb difference with: '+use_basis)
        ax3.set_ylim(-1.5,0.5)
        ax3.grid()
        ax3.legend(fontsize='xx-small')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def compare_spectra_binned_interp_mopi5_clean_factor(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', alpha=[6,7,8,9], binning=8, title='', title2='',corr_band=[]):
    fig = plt.figure(figsize=(8, 7))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].interpolated_Tb[calibration_cycle].data
    clean_f = ds_dict[use_basis].bin_freq.data
    color_spectro = {'AC240':'tab:orange', 'USRP-A':'tab:green', 'U5303':'tab:blue'}
    color_alpha = ['tab:orange','red','green','blue']
    for s in spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb[calibration_cycle].data
        #Tb_diff = (Tb-clean_Tb)/clean_Tb
        Tb_diff = Tb-clean_Tb
        ax1.plot(clean_f/1e9, ds_dict[s].interpolated_Tb[calibration_cycle].data, color=color_spectro[s], lw=0.5, label=s)
        ax1.plot(clean_f/1e9, clean_Tb, lw=0.5, color=color_spectro[use_basis], label=use_basis)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("frequency [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(2))
        ax1.set_title(title)
        ax1.grid()
        ax1.legend()
        ymax = np.nanmax(ds_dict[s].interpolated_Tb_corr[calibration_cycle].data)
        # if not np.isnan(ymax):
        #     ax2.set_ylim(ymax-4,ymax+0.5)
        #     ax2.yaxis.set_major_locator(MultipleLocator(2))
        #     ax2.yaxis.set_minor_locator(MultipleLocator(1))
        #     ax2.xaxis.set_minor_locator(MultipleLocator(5))
        #     ax2.tick_params(axis='both', which='major', labelsize=8)
        # ax2.set_xlabel('[MHz]',fontsize='x-small')
        # ax2.grid(which='both')

        if corr_band:
            freq_left=ds_dict[s].frequencies[calibration_cycle].where(
                abs(ds_dict[s].frequencies[calibration_cycle]-corr_band[s][0])<corr_band[s][1],drop=True)
            freq_right=ds_dict[s].frequencies[calibration_cycle].where(
                abs(ds_dict[s].frequencies[calibration_cycle]-corr_band[s][2])<corr_band[s][3],drop=True)
            if (s=='AC240') or (s=='U5303'):
                b = ds_dict[s].continuum_value_line_center[calibration_cycle].data - ds_dict[s].slope_indiv[calibration_cycle].data*cal_int_obj.observation_frequency
                ax1.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,
                ds_dict[s].slope_indiv[calibration_cycle].data*ds_dict[s].frequencies[calibration_cycle].data+b, color=color_spectro[s], lw=0.5)
                ax1.axvspan(freq_left[0]/1e9, freq_left[-1]/1e9 , alpha=0.25, color=color_spectro[s])
                ax1.axvspan(freq_right[0]/1e9, freq_right[-1]/1e9, alpha=0.25, color=color_spectro[s])

        for i, al in enumerate(alpha):
            lab = r'$ \alpha$ = '+str(al) + '%'
            al = al/100
            Tb_U303_cut = ds_dict['U5303'].interpolated_Tb[calibration_cycle].where(~np.isnan(Tb))
            mod_U5303 = (1-al)*Tb_U303_cut.data + al* np.nanmean(Tb_U303_cut.data)
            clean_f_smoothed = np.convolve(clean_f, np.ones((binning,))/binning, mode='full')
            smoothed_diff = np.convolve(Tb-mod_U5303, np.ones((binning,))/binning, mode='full')
            ax2.plot(clean_f_smoothed/1e9, smoothed_diff, lw=0.8, color=color_alpha[i], label=lab)
        
        ax2.set_title(title2)
        ax2.set_ylim(-1,0.5)
        ax2.yaxis.set_major_locator(MultipleLocator(0.2))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax2.set_ylabel('$\Delta T_B$ [K]')
        ax2.set_xlabel("frequency [GHz]")
        ax2.set_xlim(110.25, 111.4)
        #ax23.set_ylabel('[%]')
        ax2.grid(which='both')
        ax2.legend()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def compare_spectra_binned_interp_mopi5_clean_factor_variable(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', alpha=[6,7,8,9], binning=8, title='', title2='',broadband_bias=[]):
    fig = plt.figure(figsize=(8, 7))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].interpolated_Tb[calibration_cycle].data
    clean_f = ds_dict[use_basis].bin_freq.data
    color_alpha = ['tab:orange','red','green','blue']
    for s in spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb[calibration_cycle].data
        #Tb_diff = (Tb-clean_Tb)/clean_Tb
        Tb_diff = Tb-clean_Tb
        ax1.plot(clean_f/1e9, ds_dict[s].interpolated_Tb[calibration_cycle].data, color=color_spectro[s], lw=0.5, label=s)
        ax1.plot(clean_f/1e9, clean_Tb, lw=0.5, color=color_spectro[use_basis], label=use_basis)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("frequency [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(2))
        ax1.set_title(title)
        ax1.grid()
        
        ymax = np.nanmax(ds_dict[s].interpolated_Tb_corr[calibration_cycle].data)
        # if not np.isnan(ymax):
        #     ax2.set_ylim(ymax-4,ymax+0.5)
        #     ax2.yaxis.set_major_locator(MultipleLocator(2))
        #     ax2.yaxis.set_minor_locator(MultipleLocator(1))
        #     ax2.xaxis.set_minor_locator(MultipleLocator(5))
        #     ax2.tick_params(axis='both', which='major', labelsize=8)
        # ax2.set_xlabel('[MHz]',fontsize='x-small')
        # ax2.grid(which='both')

        Tb_U303_cut = ds_dict['U5303'].interpolated_Tb[calibration_cycle].where(~np.isnan(Tb))
        a = alpha[calibration_cycle]/100
        lab = r'$ \alpha$ = '+f'{alpha[calibration_cycle]:.0f}%'

        #Tb_corrected_factor = Tb*(1+a) - np.nanmean(Tb.data)*(a) - broadband_bias[calibration_cycle]
        
        # 
        mod_U5303 = (1-a)*Tb_U303_cut.data + a*np.nanmean(Tb_U303_cut.data)
        Tb_corrected_factor_non_lin = mod_U5303+broadband_bias[calibration_cycle]

        U5303_non_lin = Tb_U303_cut.data+broadband_bias[calibration_cycle]

        #Tb_corrected_factor = U5303_non_lin*(1-a) + a*np.nanmean(U5303_non_lin)
        #Tb_corrected_factor = Tb/(1-a) - (a/(1-a))*(np.nanmean(Tb.data) - broadband_bias[calibration_cycle])

        #Tb_diff_corrected = Tb_corrected_factor-clean_Tb
        Tb_diff_corrected = Tb-mod_U5303
        Tb_diff_corrected_non_lin_only = Tb-U5303_non_lin
        #Tb_diff_corrected_non_lin = Tb-Tb_corrected_factor
        Tb_diff_corrected_non_lin_1 = Tb-Tb_corrected_factor_non_lin

        #ax1.plot(clean_f/1e9, Tb_corrected_factor, lw=0.5, color='r', label='AC240 corrected')
        ax1.legend()

        #(1-al)*Tb_U303_cut.data + al* np.nanmean(Tb_U303_cut.data)

        clean_f_smoothed = np.convolve(clean_f, np.ones((binning,))/binning, mode='full')
        smoothed_diff_simple = np.convolve(Tb_diff, np.ones((binning,))/binning, mode='full')
        smoothed_diff = np.convolve(Tb_diff_corrected, np.ones((binning,))/binning, mode='full')
        #smoothed_diff_corr_non_lin = np.convolve(Tb_diff_corrected_non_lin, np.ones((binning,))/binning, mode='full')
        smoothed_diff_only_non_lin = np.convolve(Tb_diff_corrected_non_lin_only, np.ones((binning,))/binning, mode='full') 
        smoothed_diff_non_lin_corr = np.convolve(Tb_diff_corrected_non_lin_1, np.ones((binning,))/binning, mode='full') 
        
        
        ax2.plot(clean_f_smoothed/1e9, smoothed_diff_simple, lw=0.8, color=color_alpha[0], label=r'$ \alpha$ = 0%')
        ax2.plot(clean_f_smoothed/1e9, smoothed_diff, lw=0.8, color='r', label=lab)
        #ax2.plot(clean_f_smoothed/1e9, smoothed_diff_only_non_lin , lw=0.8, color=color_alpha[3], label=r'$ \alpha$ = 0%, with $\Delta T_{B,nonlin}$')
        
        ax2.plot(clean_f_smoothed/1e9, smoothed_diff_non_lin_corr , lw=0.8, color='g', label=r'$ \alpha$ = 8%, with $\Delta T_{B,nonlin}$')
        #ax2.plot(clean_f_smoothed/1e9, smoothed_diff_non_lin_corr , lw=0.3, color='m', label=r'$ \alpha$ = 8%,')


        ax2.axhline(0,lw=0.6, color='k', ls='--')
        ax2.set_title(title2)
        ax2.set_ylim(-1,0.4)
        ax2.yaxis.set_major_locator(MultipleLocator(0.2))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax2.set_ylabel('$\Delta T_B$ [K]')
        ax2.set_xlabel("frequency [GHz]")
        ax2.set_xlim(110.25, 111.4)
        #ax23.set_ylabel('[%]')
        ax2.grid(which='both')
        ax2.legend(loc=4)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def compare_spectra_binned_interp_mopi5_clean_factor_variable_paper(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', alpha=[6,7,8,9], binning=8, title='', title2='',broadband_bias=[]):
    fig = plt.figure(figsize=(8, 7))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].interpolated_Tb[calibration_cycle].data
    clean_f = ds_dict[use_basis].bin_freq.data
    color_alpha = ['tab:orange','red','green','blue']
    for s in spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb[calibration_cycle].data
        #Tb_OG = ds_dict[s].interpolated_Tb[calibration_cycle].data
        #Tb_diff = (Tb-clean_Tb)/clean_Tb
        Tb_diff = Tb-clean_Tb
        ax1.plot(clean_f/1e9, Tb, color=color_spectro[s], lw=0.5, label=s)
        ax1.plot(clean_f/1e9, clean_Tb, lw=0.5, color=color_spectro[use_basis], label=use_basis)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("frequency [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(2))
        ax1.set_title(title)
        ax1.grid()
        
        ymax = np.nanmax(ds_dict[s].interpolated_Tb_corr[calibration_cycle].data)
        # if not np.isnan(ymax):
        #     ax2.set_ylim(ymax-4,ymax+0.5)
        #     ax2.yaxis.set_major_locator(MultipleLocator(2))
        #     ax2.yaxis.set_minor_locator(MultipleLocator(1))
        #     ax2.xaxis.set_minor_locator(MultipleLocator(5))
        #     ax2.tick_params(axis='both', which='major', labelsize=8)
        # ax2.set_xlabel('[MHz]',fontsize='x-small')
        # ax2.grid(which='both')

        Tb_U303_cut = ds_dict['U5303'].interpolated_Tb[calibration_cycle].where(~np.isnan(Tb))
        a = alpha[calibration_cycle]/100
        lab = r'$ \alpha$ = '+f'{alpha[calibration_cycle]:.0f}%'

        #Tb_corrected_factor = Tb*(1+a) - np.nanmean(Tb.data)*(a) - broadband_bias[calibration_cycle]
        
        # 
        #mod_U5303 = (1-a)*Tb_U303_cut.data + a*np.nanmean(Tb_U303_cut.data)
        #Tb_corrected_factor_non_lin = mod_U5303+broadband_bias[calibration_cycle]

        #U5303_non_lin = Tb_U303_cut.data+broadband_bias[calibration_cycle]

        Tb_corr = (1/(1-a))*(Tb - a*np.nanmean(Tb))
        Tb_corr_all = (1/(1-a))*(Tb - a*np.nanmean(Tb) - broadband_bias[calibration_cycle])

        #Tb_corrected_factor = U5303_non_lin*(1-a) + a*np.nanmean(U5303_non_lin)
        #Tb_corrected_factor = Tb/(1-a) - (a/(1-a))*(np.nanmean(Tb.data) - broadband_bias[calibration_cycle])

        #Tb_diff_corrected = Tb_corrected_factor-clean_Tb
        #Tb_diff_corrected = Tb-mod_U5303
        #Tb_diff_corrected_non_lin_only = Tb-U5303_non_lin
        #Tb_diff_corrected_non_lin = Tb-Tb_corrected_factor
        #Tb_diff_corrected_non_lin_1 = Tb-Tb_corrected_factor_non_lin

        Tb_diff_corrected = Tb_corr - clean_Tb
        Tb_diff_corrected_non_lin = Tb_corr_all - clean_Tb


        #ax1.plot(clean_f/1e9, Tb_corrected_factor, lw=0.5, color='r', label='AC240 corrected')
        ax1.legend()

        #(1-al)*Tb_U303_cut.data + al* np.nanmean(Tb_U303_cut.data)

        clean_f_smoothed = np.convolve(clean_f, np.ones((binning,))/binning, mode='full')
        smoothed_diff_simple = np.convolve(Tb_diff, np.ones((binning,))/binning, mode='full')
        smoothed_diff = np.convolve(Tb_diff_corrected, np.ones((binning,))/binning, mode='full')
        #smoothed_diff_corr_non_lin = np.convolve(Tb_diff_corrected_non_lin, np.ones((binning,))/binning, mode='full')
        #smoothed_diff_only_non_lin = np.convolve(Tb_diff_corrected_non_lin_only, np.ones((binning,))/binning, mode='full') 
        smoothed_diff_non_lin_corr = np.convolve(Tb_diff_corrected_non_lin, np.ones((binning,))/binning, mode='full') 
        
        
        ax2.plot(clean_f_smoothed/1e9, smoothed_diff_simple, lw=0.8, color=color_alpha[0], label=r'$ \alpha$ = 0%')
        ax2.plot(clean_f_smoothed/1e9, smoothed_diff, lw=0.8, color='r', label=lab)
        #ax2.plot(clean_f_smoothed/1e9, smoothed_diff_only_non_lin , lw=0.8, color=color_alpha[3], label=r'$ \alpha$ = 0%, with $\Delta T_{B,nonlin}$')
        
        ax2.plot(clean_f_smoothed/1e9, smoothed_diff_non_lin_corr , lw=0.8, color='g', label=r'$ \alpha$ = 8%, with $\Delta T_{B,nonlin}$')
        #ax2.plot(clean_f_smoothed/1e9, smoothed_diff_non_lin_corr , lw=0.3, color='m', label=r'$ \alpha$ = 8%,')


        ax2.axhline(0,lw=0.6, color='k', ls='--')
        ax2.set_title(title2)
        ax2.set_ylim(-1,0.4)
        ax2.yaxis.set_major_locator(MultipleLocator(0.2))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax2.set_ylabel('$\Delta T_B$ [K]')
        ax2.set_xlabel("frequency [GHz]")
        ax2.set_xlim(110.25, 111.4)
        #ax23.set_ylabel('[%]')
        ax2.grid(which='both')
        ax2.legend(loc=4)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def compare_spectra_binned_interp_mopi5_clean_corr(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', title=''):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = ax1.inset_axes([0.1, 0.5, 0.2, 0.45])
    ax3 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].interpolated_Tb_corr[calibration_cycle].data
    clean_f = ds_dict[use_basis].bin_freq.data
    color_spectro = {'AC240':'tab:orange', 'USRP-A':'tab:green', 'U5303':'tab:blue'}
    for s in cal_int_obj.spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb_corr[calibration_cycle].data
        #Tb_diff = (Tb-clean_Tb)/clean_Tb
        Tb_diff = Tb-clean_Tb
        ax1.plot(clean_f/1e9, ds_dict[s].interpolated_Tb_corr[calibration_cycle].data, lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(4))
        ax1.set_title(title)
        ax1.grid()
        ax1.legend(fontsize='xx-small')
        #ax3.legend()
        ax2.plot((clean_f-cal_int_obj.observation_frequency)/1e6, ds_dict[s].interpolated_Tb_corr[calibration_cycle].data, lw=0.2, color=color_spectro[s])
        ax2.set_xlim(-10, 10)
        #ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        ymax = np.nanmax(ds_dict[s].interpolated_Tb_corr[calibration_cycle].data)
        if not np.isnan(ymax):
            ax2.set_ylim(ymax-4,ymax+0.5)
            ax2.yaxis.set_major_locator(MultipleLocator(2))
            ax2.yaxis.set_minor_locator(MultipleLocator(1))
            ax2.xaxis.set_minor_locator(MultipleLocator(5))
            ax2.tick_params(axis='both', which='major', labelsize=8)
        ax2.set_xlabel('[MHz]',fontsize='x-small')
        ax2.grid(which='both')

        if s in spectrometers:
            ax3.plot(clean_f/1e9, Tb_diff, lw=0.5, label=s, color=color_spectro[s])
        ax3.set_title('$T_b$ differences with: '+use_basis)
        ax3.set_ylim(-2,0.5)
        ax3.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax3.set_ylabel('$\Delta T_B$ [K]')
        ax3.set_xlabel("f [GHz]")
        ax3.set_xlim(110.25, 111.4)
        #ax3.set_ylabel('[%]')
        ax3.grid(which='both')
        #ax3.legend(fontsize='xx-small')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def plot_level2(ds, ac, retrieval_param, title="",figures = list()):
    '''
    Plotting function directly taken from Jonas ;)
    OG can be found in retrieval.py in MOPI retrievals
    
    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    ac : TYPE
        DESCRIPTION.
    title : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    figures : TYPE
        DESCRIPTION.

    '''
    fshift_ret = None
    if retrieval_param['retrieval_quantities'] == 'o3_h2o':
        ozone_ret, h2o_ret = ac.retrieval_quantities
    elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift':
        ozone_ret, h2o_ret, fshift_ret = ac.retrieval_quantities
        #print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
        print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
    elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit':
        ozone_ret, h2o_ret, polyfit_ret, fshift_ret = ac.retrieval_quantities
        print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
        print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
    elif retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit_sinefit':
        ozone_ret, h2o_ret, polyfit_ret, fshift_ret, sinefit_ret = ac.retrieval_quantities
        print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
        print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
        print('Sinefit : ', sinefit_ret.x)
    else:
        ozone_ret,  = ac.retrieval_quantities
    
    good_channels = ds.good_channels[retrieval_param['integration_cycle']].data == 1
    f_backend = ds.frequencies[retrieval_param['integration_cycle']].values[good_channels]
    y = ds.Tb[retrieval_param['integration_cycle']].values[good_channels]
    #y = ac.y[0]
    yf = ac.yf[0]
    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")
    
    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, y, label="observed")
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, yf, label="fitted")
    #axs[0].set_ylim(-5, 50)
    axs[0].legend()
    if fshift_ret is not None:
        axs[0].text(
            0.02,
            0.8,
            'fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3),
            transform=axs[0].transAxes,
            verticalalignment="bottom",
            horizontalalignment="left",
    )

    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
    axs[1].plot((f_backend- retrieval_param['obs_freq']) / 1e6, r_smooth, label="residuals smooth")
    axs[1].set_ylim(-0.5,0.5)
    axs[1].legend()
    axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(retrieval_param['obs_freq'] / 1e9))

    for ax in axs:
        ax.set_ylabel("$T_B$ [K]")
        ax.set_xlim([min((f_backend - retrieval_param['obs_freq']) / 1e6), max((f_backend - retrieval_param['obs_freq']) / 1e6)])
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    fig, axs = plt.subplots(1, 2, sharey=True)
    axs[0].plot(
        ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label="apriori")
    axs[0].set_xlim(-0.2,11)
    axs[0].set_xlabel("Ozone VMR [ppm]")
    axs[0].set_ylabel("Altitude [km]")
    axs[0].legend()
    axs[1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
    axs[1].set_xlabel("Measurement response")

    axs[0].grid(True)
    axs[1].grid(True)
    figures.append(fig)

    fig, axs = plt.subplots(1, 2, sharey=True)    

    axs[0].plot(ozone_ret.es * 1e6, ozone_ret.z_grid / 1e3, label="smoothing error")
    axs[0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3, label="obs error")
    axs[0].set_xlabel("$e$ [ppm]")
    axs[0].set_ylabel("Altitude [km]")
    axs[0].legend()

    # axs[1].plot(100*(ozone_ret.x - og_ozone)/og_ozone, ozone_ret.z_grid / 1e3, label="retrieval-og")
    # axs[1].plot(100*(ozone_ret.xa - og_ozone)/og_ozone, z_og / 1e3, label="apriori-og")
    # axs[0].set_xlabel("Rel diff [%]")
    # axs[0].set_ylabel("Altitude [km]")

    for avk in ozone_ret.avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1].plot(avk, ozone_ret.z_grid / 1e3)
    
    axs[1].set_xlabel("AVKM")
    axs[1].grid(True)
    axs[0].grid(True)

    fig.suptitle(title + " Ozone")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    #print('Atmospheric opacity :')
    #print(ac.ws.y_aux.value, ', mean: ', np.mean(ac.ws.y_aux.value))

    if retrieval_param['plot_opacities']:
        #opacities = ds.tropospheric_opacity[retrieval_param['integration_cycle']].values[good_channels]
        #plt.plot(f_backend,opacities,label='matlab' )
        fig, ax = plt.subplots(1, 1, sharey=True)
        ax.plot(f_backend,ac.ws.y_aux.value[0], label='ARTS' )
        ax.legend()

    if (retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift') or (retrieval_param['retrieval_quantities'] == 'o3_h2o_fshift_polyfit'):

        fig, axs = plt.subplots(2, 2, sharey=True)
        axs[0][0].semilogx(
            h2o_ret.x*1e6, h2o_ret.z_grid / 1e3, label="retrieved", marker="x"
        )
        axs[0][0].semilogx(h2o_ret.xa*1e6, h2o_ret.z_grid / 1e3, label="apriori")
        axs[0][0].set_xlabel("Water VMR [ppm]")
        axs[0][0].set_ylabel("Altitude [km]")
        axs[0][0].legend()

        axs[0][1].plot(h2o_ret.mr, h2o_ret.z_grid / 1e3)
        axs[0][1].set_xlabel("Measurement response")

        axs[1][0].semilogx(h2o_ret.es*1e6, h2o_ret.z_grid / 1e3, label="smoothing error")
        axs[1][0].semilogx(h2o_ret.eo*1e6, h2o_ret.z_grid / 1e3, label="obs error")
        axs[1][0].set_xlabel("$e$ [ppm]")
        axs[1][0].set_ylabel("Altitude [km]")
        axs[1][0].legend()

        for avk in h2o_ret.avkm:
            #if 0.8 <= np.sum(avk) <= 1.2:
            axs[1][1].plot(avk, h2o_ret.z_grid / 1e3)
        axs[1][1].set_xlabel("AVKM")


        axs[0][0].grid(True)
        axs[0][1].grid(True)
        axs[1][1].grid(True)
        axs[1][0].grid(True)

        axs[0][0].set_ylim(-0.5, 30)
        axs[0][1].set_ylim(-0.5, 30)
        axs[1][1].set_ylim(-0.5, 30)
        axs[1][0].set_ylim(-0.5, 30)

    #axs[0][0].grid(True)
    #axs[0][1].grid(True)
    #axs[1][1].grid(True)
    #axs[1][0].grid(True)

        fig.suptitle(" Water vapor retrieval (v{})".format(1))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        figures.append(fig)
    
        temp = ac.ws.t_field_raw.value.to_xarray()
        alt = ac.ws.z_field_raw.value.to_xarray()
        fig, axs = plt.subplots(1, 2)
        axs[0].plot(temp.sel(Latitude=0, Longitude=0).data, temp.Pressure)
        axs[0].invert_yaxis()
        axs[0].set_yscale('log')
        axs[0].set_xlabel('T [K]')
        axs[0].set_ylabel('$P$ [Pa]')
    
        axs[1].plot(temp.sel(Latitude=0, Longitude=0).data, alt.sel(Latitude=0, Longitude=0).data/1e3)
        #axs[1].invert_yaxis()
        #axs[1].set_yscale('log')
        axs[1].set_xlabel('T [K]')
        axs[1].set_ylabel('$Z$ [km]')
        fig.suptitle('Raw PTZ profile')
   
        figures.append(fig)   

    
    return figures

def compare_spectra_binned_interp_mopi5_clean(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', title='', corr_band=[]):
    fig = plt.figure(figsize=(8,7))
    ax1 = fig.add_subplot(211)
    ax2 = ax1.inset_axes([0.1, 0.5, 0.2, 0.45])
    ax3 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].binned_Tb[calibration_cycle].data
    clean_f = ds_dict[use_basis].bin_freq.data
    color_spectro = {'AC240':'tab:orange', 'USRP-A':'tab:green', 'U5303':'tab:blue'}
    for s in cal_int_obj.spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb[calibration_cycle].data
        #Tb_diff = (Tb-clean_Tb)/clean_Tb
        Tb_diff = Tb-clean_Tb
        ax1.plot(clean_f/1e9, ds_dict[s].interpolated_Tb[calibration_cycle].data, lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("frequency [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(4))
        ax1.set_title(title)
        ax1.grid()
        ax1.legend()
        #ax3.legend()
        ax2.plot((clean_f-cal_int_obj.observation_frequency)/1e6, ds_dict[s].interpolated_Tb[calibration_cycle].data, lw=0.2, color=color_spectro[s])
        ax2.set_xlim(-10, 10)
        #ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        ymax = np.nanmax(ds_dict[s].interpolated_Tb[calibration_cycle].data)
        if not np.isnan(ymax):
            ax2.set_ylim(ymax-4,ymax+0.5)
            ax2.yaxis.set_major_locator(MultipleLocator(2))
            ax2.yaxis.set_minor_locator(MultipleLocator(1))
            ax2.xaxis.set_minor_locator(MultipleLocator(5))
            ax2.tick_params(axis='both', which='major', labelsize=10)
        ax2.set_xlabel('[MHz]',fontsize='medium')
        ax2.grid(which='both')

        if s in spectrometers:
            ax3.plot(clean_f/1e9, Tb_diff, lw=0.5, label=s, color=color_spectro[s])
            # print('$T_b$ differences, '+s)
            # print(np.nanmean(Tb_diff))
            ax3.axhline(0,lw=0.6, color='k', ls='--')
        ax3.set_title('$T_b$ differences with: '+use_basis)
        ax3.set_ylim(-1.5,0.5)
        ax3.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax3.set_ylabel('$\Delta T_B$ [K]')
        ax3.set_xlabel("frequency [GHz]")
        ax3.set_xlim(110.25, 111.4)
        #ax3.set_ylabel('[%]')
        ax3.grid(which='both')
        #ax3.legend(fontsize='xx-small')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def compare_spectra_binned_interp_mopi5_corrected(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', title=''):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = ax1.inset_axes([0.1, 0.6, 0.15, 0.35])
    ax3 = fig.add_subplot(212)
    clean_Tb = ds_dict[use_basis].interpolated_Tb_corr[calibration_cycle].data
    clean_f = ds_dict[use_basis].bin_freq_interp.data
    for s in spectrometers:
        #mask = ds_dict[s].good_channels[calibration_cycle].data
        #mask[mask==0]=np.nan
        Tb =  ds_dict[s].interpolated_Tb_corr[calibration_cycle].data
        Tb_diff = Tb-clean_Tb
        ax1.plot(clean_f/1e9, ds_dict[s].interpolated_Tb_corr[calibration_cycle].data, lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[id].data)-10,np.median(ds_dict[s].Tb[id].data)+15)
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.set_title(title)
        ax1.grid()
        #ax1.legend(fontsize='xx-small')
        #ax3.legend()
        ax2.plot((clean_f-cal_int_obj.observation_frequency)/1e6, ds_dict[s].interpolated_Tb_corr[calibration_cycle].data, lw=0.2)
        ax2.set_xlim(-10, 10)
        #ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        ymax = np.nanmax(ds_dict[s].interpolated_Tb_corr[calibration_cycle].data)
        if not np.isnan(ymax):
            ax2.set_ylim(ymax-4,ymax+0.5)
        ax2.grid()
        ax3.plot(clean_f/1e9, Tb_diff, lw=0.5, label=s)
        ax3.set_title('Binned and interpolated Tb difference with: '+use_basis)
        ax3.set_ylim(-2,0.5)
        ax3.set_xlabel("f [GHz]")
        ax3.grid()
        ax3.legend(fontsize='xx-small')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
    return fig

def test_func(x,a,b):
    return a*np.sin(b*x)

def plot_level2_from_tropospheric_corrected_mopi5(ds, ac, retrieval_param, title="", figures=list(), fshift=True, compute_bl=True):
    '''
    Plotting function directly taken from Jonas ;)
    OG can be found in retrieval.py in MOPI retrievals
    
    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    ac : TYPE
        DESCRIPTION.
    title : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    figures : TYPE
        DESCRIPTION.

    ''' 
    from scipy import optimize
    fshift_ret = None
    compute_bl = False
    fshift = False
    if retrieval_param['retrieval_quantities'] == 'o3':
        ozone_ret, = ac.retrieval_quantities
    elif retrieval_param['retrieval_quantities'] == 'o3_fshift':
        ozone_ret, fshift_ret = ac.retrieval_quantities
        #print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
        print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
    elif retrieval_param['retrieval_quantities'] == 'o3_fshift_polyfit':
        ozone_ret, polyfit_ret, fshift_ret = ac.retrieval_quantities
        compute_bl = True
        fshift = True
        print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
        print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
    elif retrieval_param['retrieval_quantities'] == 'o3_fshift_polyfit_sinefit':
        ozone_ret, polyfit_ret, fshift_ret, sinefit_ret = ac.retrieval_quantities
        compute_bl = True
        fshift = True
        print('Poly coefficients: ' + ', '.join(['{:.2f}'.format(x[0]) for x in polyfit_ret.x]))
        print('fshift fit: {:g} kHz'.format(fshift_ret.x[0]/1e3))
        print('Sinefit : ', sinefit_ret.x)
    
    good_channels = ds.good_channels[retrieval_param['integration_cycle']].data == 1
    f_backend = ds.frequencies[retrieval_param['integration_cycle']].values[good_channels]
    y = ds.Tb_corr[retrieval_param['integration_cycle']].values[good_channels]
    #y = ac.y[0]
    yf = ac.yf[0]
    if compute_bl:
        bl = ac.y_baseline[0]

    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")
    params, params_cov = optimize.curve_fit(test_func, f_backend, r_smooth, p0=[1e-6, 2*np.pi/320e6])
    print('Fitted period on residuals: ', str(2*1e-6*np.pi/params[1]), ' MHz')

    if fshift & compute_bl:
        # Text
        fshift_text = "$\\Delta f =$ {:g} kHz".format(fshift_ret.x[0] / 1e3)
        #fshift_text='without'
        baseline_text = ", ".join(
            ["$b_{}={:.2f}$".format(i, b[0]) for i, b in enumerate(polyfit_ret.x)]
        )

    #figures = list()

    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, y, label="observed")
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, yf, label="fitted")
    if compute_bl:
        axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, bl, label="baseline")
    axs[0].set_ylim(-5, 35)
    axs[0].legend(loc='upper right')
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r_smooth, label="residuals smooth")
    axs[1].set_ylim(-1, 1)
    axs[1].legend()
    axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(retrieval_param['obs_freq'] / 1e9))
    if fshift or compute_bl:
        axs[0].text(
            0.02,
            0.8,
            fshift_text + "\n" + baseline_text,
            transform=axs[0].transAxes,
            verticalalignment="top",
            horizontalalignment="left",
        )

    for ax in axs:
        ax.set_ylabel("$T_B$ [K]")
        ax.set_xlim([min((f_backend - retrieval_param['obs_freq']) / 1e6), max((f_backend - retrieval_param['obs_freq']) / 1e6)])
    fig.suptitle(title + " - Spectrum")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    fig, axs = plt.subplots(1, 2, sharey=True)
    axs[0].plot(
        ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label="apriori")
    axs[0].set_xlim(-0.5,12)
    axs[0].set_xlabel("Ozone VMR [ppm]")
    axs[0].set_ylabel("Altitude [km]")
    axs[0].legend()
    axs[1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
    axs[1].set_xlabel("Measurement response")
    
    axs[0].grid(True)
    axs[1].grid(True)
    figures.append(fig)

    fig, axs = plt.subplots(1, 2, sharey=True)
    axs[0].plot(ozone_ret.es * 1e6, ozone_ret.z_grid / 1e3, label="smoothing error")
    axs[0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3, label="obs error")
    axs[0].set_xlabel("$e$ [ppm]")
    axs[0].set_ylabel("Altitude [km]")
    axs[0].legend()

    for avk in ozone_ret.avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1].plot(avk, ozone_ret.z_grid / 1e3)
    axs[1].set_xlabel("AVKM")
    axs[1].grid(True)
    axs[0].grid(True)


    fig.suptitle(title + " - Ozone")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    # temp = ac.ws.t_field_raw.value.to_xarray()
    # alt = ac.ws.z_field_raw.value.to_xarray()
    # fig, axs = plt.subplots(1, 2)
    # axs[0].plot(temp.sel(Latitude=0, Longitude=0).data, temp.Pressure)
    # axs[0].invert_yaxis()
    # axs[0].set_yscale('log')
    # axs[0].set_xlabel('T [K]')
    # axs[0].set_ylabel('$P$ [Pa]')

    # axs[1].plot(temp.sel(Latitude=0, Longitude=0).data, alt.sel(Latitude=0, Longitude=0).data/1e3)
    # #axs[1].invert_yaxis()
    # #axs[1].set_yscale('log')
    # axs[1].set_xlabel('T [K]')
    # axs[1].set_ylabel('$Z$ [km]')
    # fig.suptitle('Raw PTZ profile')
   
    # figures.append(fig)   

    return figures


def plot_O3_all_mopi5(level2_data, outName):
    figures2=list()
    #fig, axs = plt.subplots(nrows=3, ncols=5, sharex=True, sharey=True)
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    fig_mr = plt.figure()
    fig_mr.subplots_adjust(hspace=0.4, wspace=0.4)

    for spectro in level2_data.spectrometers:
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
            mr = level2_data[spectro].isel(observation=i, o3_lat=0, o3_lon=0).o3_m  
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
    save_single_pdf(basename_lvl2+outName+integration_strategy+'.pdf',figures2)  

def plot_O3_sel_mopi5(level2_data, spectro, outName):
    # fig = plt.figure(figsize=(9,6))
    # ax1 = fig.add_subplot(1,3,1)
    # ax2 = fig.add_subplot(1,3,2)
    # ax3 = fig.add_subplot(1,3,3   
    spectro_lvl2 = spectro
    figure_o3_sel=list()
    for i in range(len(level2_data['AC240'].observation)):
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(9,6))
        for spectro in spectro_lvl2:
            o3 = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x
            o3_apriori = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_xa
            o3_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_z
            o3_p = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_p
            mr = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
            eo = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo
            #error = lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo +  lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es
            error = np.sqrt(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es**2)
            error_frac = error/o3
            o3_good = o3.where(mr>0.8).data
            axs[0].plot(o3*1e6, o3_z/1e3, '--', linewidth=0.2, color=color_spectro[spectro])


            #axs[0].errorbar(o3*1e6, o3_z/1e3, xerr=error.values*1e6, ls='--', elinewidth=0.2, capsize=2, ecolor=color_spectro[spectro], linewidth=0.2, color=color_spectro[spectro])

            axs[0].plot(o3_good*1e6, o3_z/1e3, linewidth=1.1, label=spectro,color=color_spectro[spectro])
            axs[0].plot(o3_apriori*1e6, o3_z/1e3, '-', linewidth=0.4, label=spectro,color='k')
            axs[0].fill_betweenx(o3_z/1e3, (o3-error)*1e6,(o3+error)*1e6, color=color_spectro[spectro], alpha=0.1)

            # axs[0].errorbar(o3*1e6, o3_z/1e3, xerr=error.values*1e6, ls='--', elinewidth=0.2, capsize=2, ecolor=color_spectro[spectro], linewidth=0.2, color=color_spectro[spectro])

            # axs[0].plot(o3_good*1e6, o3_z/1e3, linewidth=1, label=spectro,color=color_spectro[spectro])
            # axs[0].plot(o3_apriori*1e6, o3_z/1e3, '--', linewidth=0.2, label=spectro,color='r')

            axs[0].set_title('$O_3$ VMR')
            axs[0].set_xlim(-0.5,11)
            axs[0].set_xlabel('$O_3$ VMR [ppm]')
            axs[0].yaxis.set_major_locator(MultipleLocator(10))
            axs[0].yaxis.set_minor_locator(MultipleLocator(5))
            axs[0].xaxis.set_major_locator(MultipleLocator(5))
            axs[0].xaxis.set_minor_locator(MultipleLocator(1))
            axs[0].grid(which='both',  axis='x', linewidth=0.5)
            axs[0].set_ylabel('Altitude [km]')

            mr_basis = level2_data['U5303'].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
            axs[1].plot(mr, o3_z/1e3, linewidth=1, color=color_spectro[spectro], label=spectro)
            axs[1].set_xlim(-0.1,1.2)
            axs[1].yaxis.set_major_locator(MultipleLocator(10))
            axs[1].yaxis.set_minor_locator(MultipleLocator(5))
            axs[1].xaxis.set_major_locator(MultipleLocator(0.4))
            axs[1].xaxis.set_minor_locator(MultipleLocator(0.1))
            axs[1].set_xlabel('MR [-]')
            axs[1].grid(which='both',  axis='x', linewidth=0.5)
            axs[1].set_title('Measurement response') 
            o3_diff = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x - level2_data['U5303'].isel(time=i, o3_lat=0, o3_lon=0).o3_x
            o3_good_diff = o3_diff.where((mr.data >0.8) & (mr_basis.data>0.8))
            if not spectro=='U5303':
                axs[2].plot(o3_diff*1e6, o3_z/1e3, '--', linewidth=0.4, color=color_spectro[spectro], label=spectro)
                axs[2].plot(o3_good_diff*1e6, o3_z/1e3, linewidth=1, color=color_spectro[spectro], label=spectro)
            axs[2].set_title('$O_3$ VMR difference with U5303')
            axs[2].set_xlim(-2,2)
            axs[2].set_xlabel('$\Delta O_3$ [ppm]')
            axs[2].yaxis.set_major_locator(MultipleLocator(10))
            axs[2].yaxis.set_minor_locator(MultipleLocator(5))
            axs[2].xaxis.set_major_locator(MultipleLocator(1))
            axs[2].xaxis.set_minor_locator(MultipleLocator(0.5))

            axs[2].grid(which='both',  axis='x', linewidth=0.5) 
            for a in axs:
                a.set_ylim(10,80)
                a.grid(which='both', axis='y', linewidth=0.5)

            # #ax2.set_xticklabels([])
            # if i==0 or i==5 or i==10:
            #     ax.set_ylabel('alt [km]', fontsize='small')
            #     ax2.set_ylabel('alt [km]', fontsize='small')
            # else: 
            #     ax.set_yticklabels([])
            #     ax2.set_yticklabels([])
            # if i in [10, 11, 12, 13, 14]:
            #     ax.set_xlabel('$O_3$ [ppm]', fontsize='small')
            #     ax2.set_xlabel('MR', fontsize='small') 
            # else: 
            #     ax.set_xticklabels([])
            #     ax2.set_xticklabels([])
        legend_elements = [
        Line2D([0], [0], color=color_spectro['U5303'], label='U5303'),
        Line2D([0], [0], color=color_spectro['AC240'], label='AC240'),
        Line2D([0], [0], color=color_spectro['USRP-A'], label='USRP-A'),
        Line2D([0], [0], color=color_spectro['AC240_unbiased'], label='AC240_unbiased'),
        Line2D([0], [0], linestyle='--', color='k', label='a priori')
        ]
        axs[0].legend(handles=legend_elements)
        axs[2].axvline(x=0, linewidth=0.6,color='k')
        fig.suptitle('$O_3$ retrievals with chunks '+str(i))
        figure_o3_sel.append(fig)
    save_single_pdf(outName+'.pdf',figure_o3_sel)


def plot_O3_chunk_mopi5(level2_data, spectro, i, title):
    # fig = plt.figure(figsize=(9,6))
    # ax1 = fig.add_subplot(1,3,1)
    # ax2 = fig.add_subplot(1,3,2)
    # ax3 = fig.add_subplot(1,3,3   
    spectro_lvl2 = spectro
    figure_o3_sel=list()
    #for i in range(len(level2_data['AC240'].observation)):
    fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(9,6))
    for spectro in spectro_lvl2:
        o3 = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x
        o3_apriori = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_xa
        o3_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_z
        o3_p = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_p
        mr = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
        eo = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo
        #error = lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo +  lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es
        error = np.sqrt(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es**2)
        error_frac = error/o3
        o3_good = o3.where(mr>0.8).data
        axs[0].plot(o3*1e6, o3_z/1e3, '--', linewidth=0.2, color=color_spectro[spectro])
        #axs[0].errorbar(o3*1e6, o3_z/1e3, xerr=error.values*1e6, ls='--', elinewidth=0.2, capsize=2, ecolor=color_spectro[spectro], linewidth=0.2, color=color_spectro[spectro])
        axs[0].plot(o3_good*1e6, o3_z/1e3, linewidth=1.1, label=spectro,color=color_spectro[spectro])
        axs[0].plot(o3_apriori*1e6, o3_z/1e3, '-', linewidth=0.4, label=spectro,color='k')
        axs[0].fill_betweenx(o3_z/1e3, (o3-error)*1e6,(o3+error)*1e6, color=color_spectro[spectro], alpha=0.1)
        axs[0].set_title('$O_3$ VMR')
        axs[0].set_xlim(-0.5,11)
        axs[0].set_xlabel('$O_3$ VMR [ppm]')
        axs[0].yaxis.set_major_locator(MultipleLocator(10))
        axs[0].yaxis.set_minor_locator(MultipleLocator(5))
        axs[0].xaxis.set_major_locator(MultipleLocator(5))
        axs[0].xaxis.set_minor_locator(MultipleLocator(1))
        axs[0].grid(which='both',  axis='x', linewidth=0.5)
        axs[0].set_ylabel('Altitude [km]')
        mr_basis = level2_data['U5303'].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
        axs[1].plot(mr, o3_z/1e3, linewidth=1, color=color_spectro[spectro], label=spectro)
        axs[1].set_xlim(-0.1,1.2)
        axs[1].yaxis.set_major_locator(MultipleLocator(10))
        axs[1].yaxis.set_minor_locator(MultipleLocator(5))
        axs[1].xaxis.set_major_locator(MultipleLocator(0.4))
        axs[1].xaxis.set_minor_locator(MultipleLocator(0.1))
        axs[1].set_xlabel('MR [-]')
        axs[1].grid(which='both',  axis='x', linewidth=0.5)
        axs[1].set_title('Measurement response') 
        o3_diff = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x - level2_data['U5303'].isel(time=i, o3_lat=0, o3_lon=0).o3_x
        o3_good_diff = o3_diff.where((mr.data >0.8) & (mr_basis.data>0.8))
        if not spectro=='U5303':
            axs[2].plot(o3_diff*1e6, o3_z/1e3, '--', linewidth=0.4, color=color_spectro[spectro], label=spectro)
            axs[2].plot(o3_good_diff*1e6, o3_z/1e3, linewidth=1, color=color_spectro[spectro], label=spectro)
        axs[2].set_title('$O_3$ VMR difference with U5303')
        axs[2].set_xlim(-2,2)
        axs[2].set_xlabel('$\Delta O_3$ [ppm]')
        axs[2].yaxis.set_major_locator(MultipleLocator(10))
        axs[2].yaxis.set_minor_locator(MultipleLocator(5))
        axs[2].xaxis.set_major_locator(MultipleLocator(1))
        axs[2].xaxis.set_minor_locator(MultipleLocator(0.5))
        axs[2].grid(which='both',  axis='x', linewidth=0.5) 
        for a in axs:
            a.set_ylim(10,80)
            a.grid(which='both', axis='y', linewidth=0.5)
        # #ax2.set_xticklabels([])
        # if i==0 or i==5 or i==10:
        #     ax.set_ylabel('alt [km]', fontsize='small')
        #     ax2.set_ylabel('alt [km]', fontsize='small')
        # else: 
        #     ax.set_yticklabels([])
        #     ax2.set_yticklabels([])
        # if i in [10, 11, 12, 13, 14]:
        #     ax.set_xlabel('$O_3$ [ppm]', fontsize='small')
        #     ax2.set_xlabel('MR', fontsize='small') 
        # else: 
        #     ax.set_xticklabels([])
        #     ax2.set_xticklabels([])
    legend_elements = [
    Line2D([0], [0], color=color_spectro['U5303'], label='U5303'),
    Line2D([0], [0], color=color_spectro['AC240'], label='AC240'),
    Line2D([0], [0], color=color_spectro['USRP-A'], label='USRP-A'),
    Line2D([0], [0], color=color_spectro['AC240_unbiased'], label='AC240_unbiased'),
    Line2D([0], [0], linestyle='--', color='k', label='a priori')
    ]
    axs[0].legend(handles=legend_elements)
    axs[2].axvline(x=0, linewidth=0.6,color='k')
    fig.suptitle(title)
    #figure_o3_sel.append(fig)
    #save_single_pdf(outName+'.pdf',figure_o3_sel)
    return fig

def plot_O3_chunk_mopi5(level2_data, spectro, i, title):
    # fig = plt.figure(figsize=(9,6))
    # ax1 = fig.add_subplot(1,3,1)
    # ax2 = fig.add_subplot(1,3,2)
    # ax3 = fig.add_subplot(1,3,3   
    spectro_lvl2 = spectro
    figure_o3_sel=list()
    #for i in range(len(level2_data['AC240'].observation)):
    fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(9,6))
    for spectro in spectro_lvl2:
        o3 = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x
        o3_apriori = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_xa
        o3_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_z
        o3_p = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_p
        mr = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
        eo = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo
        #error = lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo +  lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es
        error = np.sqrt(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es**2)
        error_frac = error/o3
        o3_good = o3.where(mr>0.8).data
        axs[0].plot(o3*1e6, o3_z/1e3, '--', linewidth=0.2, color=color_spectro[spectro])
        #axs[0].errorbar(o3*1e6, o3_z/1e3, xerr=error.values*1e6, ls='--', elinewidth=0.2, capsize=2, ecolor=color_spectro[spectro], linewidth=0.2, color=color_spectro[spectro])
        axs[0].plot(o3_good*1e6, o3_z/1e3, linewidth=1.1, label=spectro,color=color_spectro[spectro])
        axs[0].plot(o3_apriori*1e6, o3_z/1e3, '-', linewidth=0.4, label=spectro,color='k')
        axs[0].fill_betweenx(o3_z/1e3, (o3-error)*1e6,(o3+error)*1e6, color=color_spectro[spectro], alpha=0.1)
        axs[0].set_title('$O_3$ VMR')
        axs[0].set_xlim(-0.5,11)
        axs[0].set_xlabel('$O_3$ VMR [ppm]')
        axs[0].yaxis.set_major_locator(MultipleLocator(10))
        axs[0].yaxis.set_minor_locator(MultipleLocator(5))
        axs[0].xaxis.set_major_locator(MultipleLocator(5))
        axs[0].xaxis.set_minor_locator(MultipleLocator(1))
        axs[0].grid(which='both',  axis='x', linewidth=0.5)
        axs[0].set_ylabel('Altitude [km]')
        mr_basis = level2_data['U5303'].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
        axs[1].plot(mr, o3_z/1e3, linewidth=1, color=color_spectro[spectro], label=spectro)
        axs[1].set_xlim(-0.1,1.2)
        axs[1].yaxis.set_major_locator(MultipleLocator(10))
        axs[1].yaxis.set_minor_locator(MultipleLocator(5))
        axs[1].xaxis.set_major_locator(MultipleLocator(0.4))
        axs[1].xaxis.set_minor_locator(MultipleLocator(0.1))
        axs[1].set_xlabel('MR [-]')
        axs[1].grid(which='both',  axis='x', linewidth=0.5)
        axs[1].set_title('Measurement response') 
        o3_diff = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x - level2_data['U5303'].isel(time=i, o3_lat=0, o3_lon=0).o3_x
        o3_good_diff = o3_diff.where((mr.data >0.8) & (mr_basis.data>0.8))
        if not spectro=='U5303':
            axs[2].plot(o3_diff*1e6, o3_z/1e3, '--', linewidth=0.4, color=color_spectro[spectro], label=spectro)
            axs[2].plot(o3_good_diff*1e6, o3_z/1e3, linewidth=1, color=color_spectro[spectro], label=spectro)
        axs[2].set_title('$O_3$ VMR difference with U5303')
        axs[2].set_xlim(-2,2)
        axs[2].set_xlabel('$\Delta O_3$ [ppm]')
        axs[2].yaxis.set_major_locator(MultipleLocator(10))
        axs[2].yaxis.set_minor_locator(MultipleLocator(5))
        axs[2].xaxis.set_major_locator(MultipleLocator(1))
        axs[2].xaxis.set_minor_locator(MultipleLocator(0.5))
        axs[2].grid(which='both',  axis='x', linewidth=0.5) 
        for a in axs:
            a.set_ylim(10,80)
            a.grid(which='both', axis='y', linewidth=0.5)
        # #ax2.set_xticklabels([])
        # if i==0 or i==5 or i==10:
        #     ax.set_ylabel('alt [km]', fontsize='small')
        #     ax2.set_ylabel('alt [km]', fontsize='small')
        # else: 
        #     ax.set_yticklabels([])
        #     ax2.set_yticklabels([])
        # if i in [10, 11, 12, 13, 14]:
        #     ax.set_xlabel('$O_3$ [ppm]', fontsize='small')
        #     ax2.set_xlabel('MR', fontsize='small') 
        # else: 
        #     ax.set_xticklabels([])
        #     ax2.set_xticklabels([])
    legend_elements = [
    Line2D([0], [0], color=color_spectro['U5303'], label='U5303'),
    Line2D([0], [0], color=color_spectro['AC240'], label='AC240'),
    Line2D([0], [0], color=color_spectro['USRP-A'], label='USRP-A'),
    Line2D([0], [0], color=color_spectro['AC240_unbiased'], label='AC240_unbiased'),
    Line2D([0], [0], linestyle='--', color='k', label='a priori')
    ]
    axs[0].legend(handles=legend_elements)
    axs[2].axvline(x=0, linewidth=0.6,color='k')
    fig.suptitle(title)
    #figure_o3_sel.append(fig)
    #save_single_pdf(outName+'.pdf',figure_o3_sel)
    return fig

def plot_O3_all(level2_data, outName, spectro, cycles=None):
    # fig = plt.figure(figsize=(9,6))
    # ax1 = fig.add_subplot(1,3,1)
    # ax2 = fig.add_subplot(1,3,2)
    # ax3 = fig.add_subplot(1,3,3      
    
    figure_o3_sel=list()

    if cycles is None:
        cycles = np.arange(len(level2_data[spectro].time))

    for i in cycles:
        f_backend = level2_data[spectro].f.data
        y = level2_data[spectro].y[i].data
        yf = level2_data[spectro].yf[i].data
        bl = level2_data[spectro].y_baseline[i].data 
        r = y - yf
        r_smooth = np.convolve(r, np.ones(128) / 128, mode="same")

        # Text
        fshift_text = "$\\Delta f =$ {:g} kHz".format(level2_data[spectro].freq_shift_x[i].values[0] / 1e3)
        baseline_text = ", ".join(
            ["$b_{}={:.2f}$".format(a, b) for a, b in enumerate(level2_data[spectro].poly_fit_x[:,i].values)]
        )


        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(9,6))
        axs[0].plot((f_backend - F0) / 1e6, y, label="observed")
        axs[0].plot((f_backend - F0) / 1e6, yf, label="fitted")
        axs[0].text(
            0.02,
            0.8,
            fshift_text + "\n" + baseline_text,
            transform=axs[0].transAxes,
            verticalalignment="top",
            horizontalalignment="left",
        )

       # axs[0].set_xlim(-0.5,15)
        axs[0].legend()
        axs[1].plot((f_backend - F0) / 1e6, r, label="residuals")
        axs[1].plot((f_backend - F0) / 1e6, r_smooth, label="residuals smooth")
        axs[1].plot((f_backend - F0) / 1e6, bl, label="baseline")
        
       # axs[1].set_ylim(-4, 4)
        axs[1].legend()
        axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(F0 / 1e9))
        axs[1].set_ylim(-0.8,0.8)

        for ax in axs:
            ax.set_ylabel("$T_B$ [K]")
            ax.set_xlim([min((f_backend - F0) / 1e6), max((f_backend - F0) / 1e6)])
        fig.suptitle('$O_3$ retrievals for '+spectro+ ' chunk: '+str(i))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        figure_o3_sel.append(fig)

        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(9,6))
        
        o3 = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x
        o3_apriori = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_xa
        o3_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_z
        o3_p = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_p
        mr = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
        #error = lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo +  lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es
        error = np.sqrt(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es**2)
        error_frac = error/o3
        o3_good = o3.where(mr>0.8).data
        #axs[0].plot(o3_good*1e6, o3_z/1e3, '--', linewidth=1, color='tab:blue')
        axs[0].plot(o3*1e6, o3_z/1e3,'-x', linewidth=1, label='retrieved',color='blue')
        axs[0].plot(o3_apriori*1e6, o3_z/1e3, '-', linewidth=0.8, label='apriori',color='r')
        axs[0].set_title('$O_3$ VMR')
        axs[0].set_xlim(-0.5,11)
        axs[0].set_ylim(min(o3_z/1e3),max(o3_z/1e3))
        axs[0].set_xlabel('$O_3$ VMR [ppm]')
        axs[0].yaxis.set_major_locator(MultipleLocator(10))
        axs[0].yaxis.set_minor_locator(MultipleLocator(5))
        axs[0].xaxis.set_major_locator(MultipleLocator(5))
        axs[0].xaxis.set_minor_locator(MultipleLocator(1))
        axs[0].grid(which='both',  axis='x', linewidth=0.5)
        axs[0].set_ylabel('Altitude [km]')
        axs[0].legend()
        axs[1].plot(mr/2, o3_z/1e3,color='k', label='MR/2')
        counter=0
        for avk in level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_avkm:
            if 0.8 <= np.sum(avk) <= 1.2:
                counter=counter+1
                if np.mod(counter,5)==0:
                    axs[1].plot(avk, o3_z / 1e3, label='z='+f'{o3_z.sel(o3_p=avk.o3_p).values/1e3:.0f}'+'km', color='r')
                else:
                    axs[1].plot(avk, o3_z / 1e3, color='k')
        axs[1].set_xlabel("AVKM")
        axs[1].set_xlim(-0.05,0.6)
        axs[1].xaxis.set_major_locator(MultipleLocator(0.1))
        axs[1].xaxis.set_minor_locator(MultipleLocator(0.05))
        axs[1].legend()
        axs[1].grid(which='both',  axis='x', linewidth=0.5)
        

        axs[2].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es * 1e6, o3_z / 1e3, label="smoothing error")
        axs[2].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo * 1e6, o3_z / 1e3, label="obs error")
        axs[2].set_xlabel("$e$ [ppm]")
        axs[2].set_ylabel("Altitude [km]")
        axs[2].legend()
        axs[2].grid(axis='x', linewidth=0.5)
        
        #axs[3].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).h2o_x * 1e6, o3_z / 1e3, label="retrieved")
        # axs[3].set_xlabel("$VMR$ [ppm]")
        # axs[3].set_ylabel("Altitude [km]")
        # axs[3].legend()
        #axs[3].grid(axis='x', linewidth=0.5)

        for a in axs:
            #a.set_ylim(10,80)
            a.grid(which='both', axis='y', linewidth=0.5)
        fig.suptitle('$O_3$ retrievals for '+spectro+ ' chunk: '+str(i))
        figure_o3_sel.append(fig)

        # #if retrieval_param['retrieval_quantities'] == 'o3_h2o':
        # fig, axs = plt.subplots(1, 1, sharey=True)
        # h2o_x = level2_data[spectro].isel(time=i).h2o-pwr98__h2o_x
        # h2o_xa = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).h2o-pwr98__h2o_xa
        # h2o_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).h2o-pwr98__h2o_z

        # axs[0].semilogx(
        #     h2o_x, h2o_z / 1e3, label="retrieved", marker="x"
        # )
        # axs[0].semilogx(h2o_xa, h2o_z / 1e3, label="apriori")
        # axs[0].set_xlabel("Water VMR []")
        # axs[0].set_ylabel("Altitude [km]")
        # axs[0].legend()
        
        # fig.suptitle(r'$H_{2}O$ retrievals (and h2o)')
        # figure_o3_sel.append(fig)
    save_single_pdf(outName+'.pdf',figure_o3_sel)

def plot_O3_sel_paper(level2_data, outName, spectro, cycles=None):
    # fig = plt.figure(figsize=(9,6))
    # ax1 = fig.add_subplot(1,3,1)
    # ax2 = fig.add_subplot(1,3,2)
    # ax3 = fig.add_subplot(1,3,3      
    
    figure_o3_sel=list()

    if cycles is None:
        cycles = np.arange(len(level2_data[spectro].time))

    for i in cycles:
        f_backend = level2_data[spectro].f.data
        y = level2_data[spectro].y[i].data
        yf = level2_data[spectro].yf[i].data
        bl = level2_data[spectro].y_baseline[i].data 
        r = y - yf
        r_smooth = np.convolve(r, np.ones(128) / 128, mode="same")

        # Text
        fshift_text = "$\\Delta f =$ {:g} kHz".format(level2_data[spectro].freq_shift_x[i].values[0] / 1e3)
        baseline_text = ", ".join(
            ["$b_{}={:.2f}$".format(a, b) for a, b in enumerate(level2_data[spectro].poly_fit_x[:,i].values)]
        )


        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(9,6))
        axs[0].plot((f_backend - F0) / 1e6, y, label="observed")
        axs[0].plot((f_backend - F0) / 1e6, yf, label="fitted")
        axs[0].text(
            0.02,
            0.8,
            fshift_text + "\n" + baseline_text,
            transform=axs[0].transAxes,
            verticalalignment="top",
            horizontalalignment="left",
        )

       # axs[0].set_xlim(-0.5,15)
        axs[0].legend()
        axs[1].plot((f_backend - F0) / 1e6, r, label="residuals")
        axs[1].plot((f_backend - F0) / 1e6, r_smooth, label="residuals smooth")
        axs[1].plot((f_backend - F0) / 1e6, bl, label="baseline")
        
       # axs[1].set_ylim(-4, 4)
        axs[1].legend()
        axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(F0 / 1e9))
        axs[1].set_ylim(-0.8,0.8)

        for ax in axs:
            ax.set_ylabel("$T_B$ [K]")
            ax.set_xlim([min((f_backend - F0) / 1e6), max((f_backend - F0) / 1e6)])
        fig.suptitle('$O_3$ retrievals for '+spectro+ ' chunk: '+str(i))
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        figure_o3_sel.append(fig)

        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(9,6))
        
        o3 = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x
        o3_apriori = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_xa
        o3_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_z
        o3_p = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_p
        mr = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
        #error = lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo +  lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es
        error = np.sqrt(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es**2)
        error_frac = error/o3
        o3_good = o3.where(mr>0.8).data
        #axs[0].plot(o3_good*1e6, o3_z/1e3, '--', linewidth=1, color='tab:blue')
       # axs[0].plot(o3*1e6, o3_z/1e3,'-x', linewidth=1, label='retrieved',color='blue')
        axs[0].plot(o3_apriori*1e6, o3_z/1e3, '-', linewidth=0.8, label='apriori',color='k')       
        axs[0].plot(o3*1e6, o3_z/1e3,'-', linewidth=1, label='retrieved',color=color_spectro[spectro])
        axs[0].fill_betweenx(o3_z/1e3, (o3-error)*1e6,(o3+error)*1e6, color=color_spectro[spectro], alpha=0.5)

       # axs[0].set_title('$O_3$ VMR')
        axs[0].set_xlim(-0.5,11)
        axs[0].set_ylim(min(o3_z/1e3),max(o3_z/1e3))
        axs[0].set_xlabel('Ozone VMR [ppmv]')
        axs[0].yaxis.set_major_locator(MultipleLocator(10))
        axs[0].yaxis.set_minor_locator(MultipleLocator(5))
        axs[0].xaxis.set_major_locator(MultipleLocator(5))
        axs[0].xaxis.set_minor_locator(MultipleLocator(1))
        axs[0].grid(which='both',  axis='x', linewidth=0.5)
        axs[0].set_ylabel('Altitude [km]')
        axs[0].legend()
        axs[1].plot(mr/2, o3_z/1e3,color='b', label='MR/2')
        counter=0
        for avk in level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_avkm:
            if 0.6 <= np.sum(avk) <= 1.4:
                counter=counter+1
                if np.mod(counter,5)==0:
                    axs[1].plot(avk, o3_z / 1e3, label='z ='+f'{o3_z.sel(o3_p=avk.o3_p).values/1e3:.0f}'+' km', color='r')
                else:
                    axs[1].plot(avk, o3_z / 1e3, color='k')
        axs[1].set_xlabel("AVK")
        axs[1].set_xlim(-0.05,0.6)
        axs[1].xaxis.set_major_locator(MultipleLocator(0.1))
        axs[1].xaxis.set_minor_locator(MultipleLocator(0.05))
        axs[1].legend()
        axs[1].grid(which='both',  axis='x', linewidth=0.5)
        

        axs[2].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es * 1e6, o3_z / 1e3, color='g', label="smoothing")
        axs[2].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo * 1e6, o3_z / 1e3, color='r', label="observation")
        axs[2].set_xlabel("Error [ppmv]")
        axs[2].set_ylabel("Altitude [km]")
        axs[2].legend(loc='upper center')
        axs[2].grid(axis='x', linewidth=0.5)
        
        #axs[3].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).h2o_x * 1e6, o3_z / 1e3, label="retrieved")
        # axs[3].set_xlabel("$VMR$ [ppm]")
        # axs[3].set_ylabel("Altitude [km]")
        # axs[3].legend()
        #axs[3].grid(axis='x', linewidth=0.5)

        for a in axs:
            #a.set_ylim(10,80)
            a.grid(which='both', axis='y', linewidth=0.5)
        fig.suptitle('Ozone retrievals for '+spectro+ ' chunk: '+str(i))
        figure_o3_sel.append(fig)

        # #if retrieval_param['retrieval_quantities'] == 'o3_h2o':
        # fig, axs = plt.subplots(1, 1, sharey=True)
        # h2o_x = level2_data[spectro].isel(time=i).h2o-pwr98__h2o_x
        # h2o_xa = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).h2o-pwr98__h2o_xa
        # h2o_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).h2o-pwr98__h2o_z

        # axs[0].semilogx(
        #     h2o_x, h2o_z / 1e3, label="retrieved", marker="x"
        # )
        # axs[0].semilogx(h2o_xa, h2o_z / 1e3, label="apriori")
        # axs[0].set_xlabel("Water VMR []")
        # axs[0].set_ylabel("Altitude [km]")
        # axs[0].legend()
        
        # fig.suptitle(r'$H_{2}O$ retrievals (and h2o)')
        # figure_o3_sel.append(fig)
    save_single_pdf(outName+'.pdf',figure_o3_sel)

def plot_O3_3on1_paper(level2_data, outName, spectrometer, cycles=[0]):
    # fig = plt.figure(figsize=(9,6))
    # ax1 = fig.add_subplot(1,3,1)
    # ax2 = fig.add_subplot(1,3,2)
    # ax3 = fig.add_subplot(1,3,3      
    
    figure_o3_sel=list()
    
    # plt.rcParams['pdf.fonttype'] = 42
    # plt.rcParams['ps.fonttype'] = 42
    #plt.rcParams['text.usetex'] = True

    for i in cycles:
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(9,6))

        for spectro in spectrometer:
            o3 = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_x
            o3_apriori = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_xa
            o3_z = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_z
            o3_p = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_p
            mr = level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_mr
            #error = lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo +  lvl2[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es
            error = np.sqrt(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo**2 +  level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es**2)
            error_frac = error/o3
            o3_good = o3.where(mr>0.8).data
            #axs[0].plot(o3_good*1e6, o3_z/1e3, '--', linewidth=1, color='tab:blue')
       #     axs[0].plot(o3*1e6, o3_z/1e3,'-x', linewidth=1, label='retrieved',color='blue')
            axs[0].plot(o3_apriori*1e6, o3_z/1e3, '-', linewidth=0.8, label='apriori',color='k')       
            axs[0].plot(o3*1e6, o3_z/1e3,'-', linewidth=1, label='retrieved',color=color_spectro[spectro])
            axs[0].fill_betweenx(o3_z/1e3, (o3-error)*1e6,(o3+error)*1e6, color=color_spectro[spectro], alpha=0.5)

       #     axs[0].set_title('$O_3$ VMR')
            axs[0].set_xlim(-0.5,11)
            axs[0].set_ylim(min(o3_z/1e3),max(o3_z/1e3))
            axs[0].set_xlabel('Ozone VMR [ppmv]')
            axs[0].yaxis.set_major_locator(MultipleLocator(10))
            axs[0].yaxis.set_minor_locator(MultipleLocator(5))
            axs[0].xaxis.set_major_locator(MultipleLocator(5))
            axs[0].xaxis.set_minor_locator(MultipleLocator(1))
            axs[0].grid(which='both',  axis='x', linewidth=0.5)
            axs[0].set_ylabel('Altitude [km]')

            legend_elements = [
                Line2D([0], [0], color='k', label='apriori'),
                Line2D([0], [0], color=color_spectro['U5303'], label='U5303'),
                Line2D([0], [0], color=color_spectro['AC240'], label='AC240'),
                Line2D([0], [0], color=color_spectro['USRP-A'], label='USRP-A'),
            ]

            #axs[0].legend()
            axs[0].legend(handles=legend_elements, fontsize=12, loc='upper right')
            
            axs[1].plot(mr, o3_z/1e3, color=color_spectro[spectro], label=spectro)
            axs[1].set_xlabel("Measurement response")
            axs[1].set_ylabel("Altitude [km]")
            axs[1].set_xlim(-0.05,1.2)
            axs[1].xaxis.set_major_locator(MultipleLocator(0.2))
            axs[1].xaxis.set_minor_locator(MultipleLocator(0.1))
            #axs[1].legend([''])
            axs[1].grid(which='both',  axis='x', linewidth=0.5)
        

            axs[2].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_es * 1e6, o3_z / 1e3, '--', color=color_spectro[spectro], label="smoothing")
            axs[2].plot(level2_data[spectro].isel(time=i, o3_lat=0, o3_lon=0).o3_eo * 1e6, o3_z / 1e3, '-', color=color_spectro[spectro], label="observation")
            axs[2].xaxis.set_major_locator(MultipleLocator(0.2))
            axs[2].xaxis.set_minor_locator(MultipleLocator(0.1))
            axs[2].set_xlabel("Errors [ppmv]")
            axs[2].set_ylabel("Altitude [km]")
            axs[2].grid(which='both',  axis='x', linewidth=0.5)


            legend_elements2 = [
                Line2D([0], [0], color='k', ls='--', label='smoothing'),
                Line2D([0], [0], color='k', ls='-',label='observation'),
            ]
            #axs[2].legend(loc='upper center')
            axs[2].legend(handles=legend_elements2, fontsize=12, loc='upper right')


            for a in axs:
                #a.set_ylim(10,80)
                a.grid(which='both', axis='y', linewidth=0.5)
        fig.suptitle('Ozone retrievals for ' + ' chunk: '+str(i))
        figure_o3_sel.append(fig)

    save_single_pdf(outName+'.pdf',figure_o3_sel)