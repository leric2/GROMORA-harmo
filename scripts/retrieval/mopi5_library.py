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

def return_bad_channels_mopi5(number_of_channel, date, spectro):
    '''
    common function for the 2 classes... maybe a better way to do it.
    '''
    if spectro == 'USRP-A':
        bad_channels = np.hstack((np.arange(0,1023), 8092, np.arange(number_of_channel-1023,number_of_channel)))
    elif spectro == 'U5303':
        #bad_channels = np.where(intermediate_frequency > 1000)
        bad_channels = np.hstack((np.arange(0,63), np.arange(11000,number_of_channel)))
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
    ax2 = ax1.inset_axes([0.1, 0.5, 0.3, 0.45])
    for s in spectrometers:
        mask = ds_dict[s].good_channels[calibration_cycle].data
        mask[mask==0]=np.nan
        ax1.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,
                 ds_dict[s].Tb[calibration_cycle].data*mask, lw=0.5, label=s)
        ax1.set_xlim(110.25, 111.4)
        #ax1.set_ylim(np.median(ds_dict[s].Tb[calibration_cycle].data)-10,np.median(ds_dict[s].Tb[id].data)+calibration_cycle)
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(2))
        ax1.set_title(title)
        
        ax1.legend(fontsize='small',loc=1)
        if corr_band:
            color_spectro = {'AC240':'tab:orange', 'USRP-A':'tab:green', 'U5303':'tab:blue'}
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
        

        #ax3.legend()
        #ax2.plot(ds_dict[s].frequencies[calibration_cycle].data/1e9,ds_dict[s].Tb[calibration_cycle].data*mask, lw=0.2)
        ax2.plot(
            (ds_dict[s].frequencies[calibration_cycle].data-cal_int_obj.observation_frequency)/1e6,
            ds_dict[s].Tb[calibration_cycle].data*mask, lw=0.2
        )
        ax2.yaxis.set_major_locator(MultipleLocator(2))
        ax2.yaxis.set_minor_locator(MultipleLocator(1))
        ax2.set_xlim(-10, 10)
        ax2.set_xlabel('[MHz]')
        #ax2.set_yticklabels([])
        #ax2.set_xticklabels([])
        ymax = np.nanmax(ds_dict[s].Tb[calibration_cycle].data*mask)
        if not np.isnan(ymax):
            ax2.set_ylim(ymax-3,ymax+1)
        
        #ax2.legend(fontsize='xx-small')
    ax1.grid()
    ax2.grid(which='both')
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

    ax32.set_ylabel('$P_{air}$ [HPa]', color='k')
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

def plot_ts_mopi5_Feb(calibration, title):
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
        dates = [pd.to_datetime(d) for d in calibration.calibrated_data[spectro].time.data]
        ax[0].scatter(dates, calibration.calibrated_data[spectro].mean_Tb, marker=marker[a], s=2, label=spectro)
        #ax[1].plot(calibration.calibrated_data[s].time, calibration.calibrated_data[s].TSys,'.',markersize=4,label=s)
        #ax2.plot(calibration.calibration_flags[s].time, calibration.calibration_flags[s].sum(dim='flags').to_array()[0,:],'.',markersize=4,label=s)
        #ax.legend()
    ax[0].set_ylim(70,260)
    ax[0].legend(fontsize='xx-small')
    s = 'AC240'

    ts1 = calibration.meteo_data[s].where(calibration.meteo_data[s].time < pd.to_datetime('2019-02-02'), drop=True)
   
    ts2 = calibration.meteo_data[s].where(calibration.meteo_data[s].time > pd.to_datetime('2019-02-02'), drop=True)
    ts2 = ts2.where(calibration.meteo_data[s].time < pd.to_datetime('2019-02-11'), drop=True)
    
    ts3 = calibration.meteo_data[s].where(calibration.meteo_data[s].time > pd.to_datetime('2019-02-11'), drop=True)

    ax31.plot(ts1.time, ts1.air_temperature-273.15,'r-', lw=meteo_lw) 
    ax31.plot(ts2.time, ts2.air_temperature-273.15,'r-', lw=meteo_lw) 
    ax31.plot(ts3.time, ts3.air_temperature-273.15,'r-', lw=meteo_lw) 

    ax32.plot(ts1.time, ts1.air_pressure,'k-', lw=meteo_lw)  
    ax32.plot(ts2.time, ts2.air_pressure,'k-', lw=meteo_lw)  
    ax32.plot(ts3.time, ts3.air_pressure,'k-', lw=meteo_lw)  
    ax31.set_ylabel('$T_{air}$ [$\degree C$]', color='r')
    ax31.set_ylim(-6, 22)
    ax31.tick_params(axis='y', labelcolor='r')

    ax32.set_ylabel('$P_{air}$ [HPa]', color='k')
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

    midmonth = mdates.DayLocator(15) 
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
    ax[1].xaxis.set_major_locator(mdates.WeekdayLocator(interval=1))
    
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

def compare_spectra_binned_interp_mopi5_clean_factor(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', alpha=[6,7,8,9], binning=8, title='',corr_band=[]):
    fig = plt.figure(figsize=(8, 11))
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
        ax1.set_xlabel("f [GHz]")
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
        
        ax2.set_title('$\Delta T_B$ using $T_B^{\'}= (1-\\alpha) T_B + \\alpha T_{B,mean}$ ')
        ax2.set_ylim(-1,0.5)
        ax2.yaxis.set_major_locator(MultipleLocator(0.2))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax2.set_ylabel('$\Delta T_B$ [K]')
        ax2.set_xlabel("f [GHz]")
        ax2.set_xlim(110.25, 111.4)
        #ax23.set_ylabel('[%]')
        ax2.grid(which='both')
        ax2.legend()
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

def compare_spectra_binned_interp_mopi5_clean(cal_int_obj, ds_dict, calibration_cycle=0, spectrometers=['AC240','USRP-A'], use_basis='U5303', title='', corr_band=[]):
    fig = plt.figure(figsize=(9,6))
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
        ax1.set_xlabel("f [GHz]")
        ax1.set_ylabel(r"$T_B$ [K]")
        ax1.yaxis.set_major_locator(MultipleLocator(4))
        ax1.set_title(title)
        ax1.grid()
        ax1.legend(fontsize='xx-small')
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
            ax2.tick_params(axis='both', which='major', labelsize=8)
        ax2.set_xlabel('[MHz]',fontsize='x-small')
        ax2.grid(which='both')

        if s in spectrometers:
            ax3.plot(clean_f/1e9, Tb_diff, lw=0.5, label=s, color=color_spectro[s])
        ax3.set_title('$T_b$ differences with: '+use_basis)
        ax3.set_ylim(-1.5,0.5)
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
    if fshift or compute_bl:
        ozone_ret, fshift_ret, polyfit_ret = ac.retrieval_quantities
    else:
        [ozone_ret] = ac.retrieval_quantities
    good_channels = ds.good_channels[retrieval_param['integration_cycle']].data == 1
    f_backend = ds.frequencies[retrieval_param['integration_cycle']].values[good_channels]
    y = ds.Tb_corr[retrieval_param['integration_cycle']].values[good_channels]
    #y = ac.y[0]
    yf = ac.yf[0]
    if compute_bl:
        bl = ac.y_baseline[0]

    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")

    if fshift:
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
    axs[0].set_ylim(-5, 50)
    axs[0].legend(loc='upper right')
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r_smooth, label="residuals smooth")
    axs[1].set_ylim(-2, 2)
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

    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].plot(
        ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlim(-0.5,15)
    axs[0][0].set_xlabel("Ozone VMR [ppm]")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")
    
    axs[1][0].plot(ozone_ret.es * 1e6, ozone_ret.z_grid / 1e3, label="smoothing error")
    axs[1][0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3, label="obs error")
    axs[1][0].set_xlabel("$e$ [ppm]")
    axs[1][0].set_ylabel("Altitude [km]")
    axs[1][0].legend()

    for avk in ozone_ret.avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1][1].plot(avk, ozone_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")

    axs[0][0].grid(True)
    axs[0][1].grid(True)
    axs[1][1].grid(True)
    axs[1][0].grid(True)


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