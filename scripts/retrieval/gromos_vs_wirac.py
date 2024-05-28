#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

A short script to compare the WIRA-C to the GROMOS level 1 spectra during the measurement campaign done at IAP.

"""
from curses import window
import sys

import datetime
import os
import time
from abc import ABC
import matplotlib

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv

from gromora_utils import save_single_pdf
load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')
import gromos_classes as gc
import wirac_classes as wc

sys.path.insert(0, '/home/esauvageat/Documents/GROMORA/Analysis/GROMORA-harmo/scripts/retrieval/')
sys.path.insert(0, '/home/esauvageat/Documents/GROMORA/Analysis/GROMORA-harmo/scripts/pyretrievals/')

from retrievals import data

ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

%matplotlib widget
OBSERVATION_FREQ = 1.4217504e11

def compare_ds(gromos, wira0, wira1, wincorr=True, time_ind=[0,1],avgN=8):
    avgN2 = avgN*4
    figures = []
    fs = 22
    for i in time_ind:
        fig, axs = plt.subplots(2, 1, figsize=(16,12))
        datstr = pd.to_datetime(gromos.isel(time=i).time.data).strftime('%Y.%m.%d:%HH')
        if wincorr:
            axs[0].plot(np.convolve(gromos.isel(time=i).frequencies/1e9, np.ones(avgN)/avgN, mode='valid') , np.convolve(gromos.isel(time=i).Tb_win_corr, np.ones(avgN)/avgN, mode='valid') , color='#1b9e77', label='GROMOS')
        else:
            axs[0].plot(np.convolve(gromos.isel(time=i).frequencies/1e9, np.ones(avgN)/avgN, mode='valid') , np.convolve(gromos.isel(time=i).Tb, np.ones(avgN)/avgN, mode='valid') , color='#1b9e77', label='GROMOS')

        axs[0].plot(np.convolve(wira0.isel(time=i).frequencies/1e9, np.ones(avgN)/avgN, mode='valid'), np.convolve(wira0.isel(time=i).Tb, np.ones(avgN)/avgN, mode='valid') , color='#d95f02', label='USRP0')
        axs[0].plot(np.convolve(wira1.isel(time=i).frequencies/1e9, np.ones(avgN)/avgN, mode='valid'), np.convolve(wira1.isel(time=i).Tb, np.ones(avgN)/avgN, mode='valid') , color='#7570b3', label='USRP1')
        
        ax0_ins = axs[0].inset_axes([0.1, 0.6, 0.15, 0.35])
        if wincorr:
            ax0_ins.plot((gromos.isel(time=i).frequencies - OBSERVATION_FREQ)/1e6 ,gromos.isel(time=i).Tb_win_corr , color='#1b9e77')
        else:
            ax0_ins.plot((gromos.isel(time=i).frequencies - OBSERVATION_FREQ)/1e6 ,gromos.isel(time=i).Tb, color='#1b9e77')

        ax0_ins.plot((wira0.isel(time=i).frequencies- OBSERVATION_FREQ )/1e6,wira0.isel(time=i).Tb_win_corr , color='#d95f02')
        ax0_ins.plot((wira1.isel(time=i).frequencies - OBSERVATION_FREQ)/1e6,wira1.isel(time=i).Tb_win_corr , color='#7570b3')
        ax0_ins.set_ylim(gromos.isel(time=i).Tb_win_corr.quantile(0.98), gromos.isel(time=i).Tb_win_corr.quantile(0.999)+5)
        ax0_ins.set_xlim(-3, 3)
        ax0_ins.set_xlabel('[MHz]')
        w0 = wira0.isel(time=i).interp_like(gromos.isel(time=i).frequencies)
        clean_Tb0 = wira0.isel(time=i).Tb
        clean_f0 = wira0.isel(time=i).frequencies

        if wincorr:
            clean_Tb_basis = gromos.isel(time=i).Tb_win_corr
        else:
            clean_Tb_basis = gromos.isel(time=i).Tb
        grid = gromos.isel(time=i).frequencies

        interp_tb0 = data.interpolate(
            grid.data, 
            clean_f0.data, 
            clean_Tb0.data, 
            left=np.nan, 
            right=np.nan
            )
            #da_interp_Tb[i,:] = interp_tb
        
        clean_Tb1 = wira1.isel(time=i).Tb
        clean_f1 = wira1.isel(time=i).frequencies

        grid = gromos.isel(time=i).frequencies

        interp_tb1 = data.interpolate(
            grid.data, 
            clean_f1.data, 
            clean_Tb1.data, 
            left=np.nan, 
            right=np.nan
            )
            #da_interp_Tb[i,:] = interp_tb
        
        #axs[0].plot(np.convolve(gromos.isel(time=i).frequencies/1e9, np.ones(avgN2)/avgN2, mode='valid'), np.convolve(interp_tb0, np.ones(avgN2)/avgN2, mode='valid') , color='r', label='USRP0 interp')
        #axs[0].plot(np.convolve(gromos.isel(time=i).frequencies/1e9, np.ones(avgN2)/avgN2, mode='valid'), np.convolve(interp_tb1, np.ones(avgN2)/avgN2, mode='valid') , color='k', label='USRP1 interp')
        axs[0].legend()
        axs[0].set_ylim(gromos.isel(time=i).Tb.median()-30,gromos.isel(time=i).Tb.median()+25)
        axs[0].set_xlim(141.7,142.7)
        axs[0].set_xlabel('RF [GHz]', fontsize=fs)
        axs[0].set_ylabel('TB [K]', fontsize=fs)
        axs[0].grid()
        axs[0].set_title('L1, '+datstr+', 1h integration', fontsize=fs)

        Tb_diff0 = clean_Tb_basis-interp_tb0
        Tb_diff1 = clean_Tb_basis-interp_tb1

        #self.integrated_data[s] = self.integrated_data[s].assign(interpolated_Tb = da_interp_Tb)

        #w1 = wira1.isel(time=i).interp_like(gromos.isel(time=i).frequencies)
        axs[1].plot(np.convolve(gromos.isel(time=i).frequencies/1e9, np.ones(avgN2)/avgN2, mode='valid'), np.convolve(Tb_diff0, np.ones(avgN2)/avgN2, mode='valid') , color='#d95f02', label='USRP0')
        axs[1].plot(np.convolve(gromos.isel(time=i).frequencies/1e9, np.ones(avgN2)/avgN2, mode='valid'), np.convolve(Tb_diff1, np.ones(avgN2)/avgN2, mode='valid') , color='#7570b3', label='USRP1')
        axs[1].legend()
        axs[1].set_ylim(-10, 10)
        axs[1].set_xlim(141.7,142.7)
        axs[1].set_xlabel('RF [GHz]', fontsize=fs)
        axs[1].set_ylabel(r'$\Delta$TB [K]', fontsize=fs)
        axs[1].grid()
        axs[1].set_title('GROMOS - WIRAC', fontsize=fs)
        plt.tight_layout(rect=[0, 0.01, 0.99, 1])
        figures.append(fig)
    return figures

if __name__ == "__main__":
    instrument_name = "GROMOS"
    instrument_name2 = "WIRAC"
    date = datetime.date(2022, 3, 2)
    int_time = 1
    integration_strategy = 'classic'
    recheck_channels = False

    windows_corrected = True

    basename_lvl2 = "/home/esauvageat/Documents/WIRAC/Level1/"
    basename_lvl1 = os.path.join(
        '/storage/atmosphere/instruments/gromos/level1/GROMORA/v2/', str(date.year))
    instrument = gc.GROMOS_LvL2(
        date,
        basename_lvl1,
        basename_lvl2,
        integration_strategy,
        int_time,
        extra_base=''
        )

    ds_gromos = instrument.read_level1b(extra_base='')[0]

    basename_lvl1_wirac = os.path.join(
        '/home/esauvageat/Documents/WIRAC/Level1/')
    wirac = wc.WIRAC_LvL2(
        date,
        basename_lvl1_wirac,
        basename_lvl2,
        integration_strategy,
        int_time,
        extra_base=''
        )
    
    ds_wirac = wirac.read_level1b(extra_base='')[0]
    ds_wirac['AC240'].Tb[:,wirac.return_bad_channels(date, 'AC240')] = np.nan
    ds_wirac['USRP0'] = ds_wirac['AC240'].where(ds_wirac['AC240'].channel_idx<=16384, drop=True)
    ds_wirac['USRP1'] = ds_wirac['AC240'].where(ds_wirac['AC240'].channel_idx>16384, drop=True)
    # instrument.plot_level1b_TB()
    # wirac.plot_level1b_TB()

    #figures = compare_ds(ds_gromos['AC240'], ds_wirac['USRP1'],ds_wirac['USRP0'], wincorr=windows_corrected, time_ind=[12,13,14,15,16,16,18,19,20,21,22,23],avgN=8)
    figures = compare_ds(ds_gromos['AC240'], ds_wirac['USRP1'],ds_wirac['USRP0'], time_ind=[5],avgN=8)
    #figures = compare_ds(ds_gromos['AC240'], ds_wirac['USRP1'],ds_wirac['USRP0'], wincorr=windows_corrected, time_ind=[0,1,2,3,4,5,6,7,8,9,10,11], avgN=8)

    # if windows_corrected:
    #     save_single_pdf(basename_lvl2+'comparison_wirac_gromos_corr_'+date.strftime('%Y.%m.%d')+'.pdf', figures)
    # else:
    #     save_single_pdf(basename_lvl2+'comparison_wirac_gromos_og_'+date.strftime('%Y.%m.%d')+'.pdf', figures)
