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

if __name__ == "__main__":
    instrument_name = "mopi5"
    #date = datetime.date(2019,2,1)
    date = pd.date_range(start='2019-01-03', end='2019-01-05')
    #date = pd.date_range(start='2019-01-30', end='2019-06-18')
    #date = pd.date_range(start='2019-01-30', end='2019-02-22')

    #date = pd.date_range(start='2019-05-01', end='2019-05-04')
    #date = pd.date_range(start='2019-04-25', end='2019-04-27')
    #date = pd.date_range(start='2019-06-11', end='2019-06-11')
    #date = pd.date_range(start='2019-03-12', end='2019-03-12')
    # options are: 'TOD', 'TOD_harmo', 'classic' 'meanTb_harmo', or 'meanTb'
    integration_strategy = 'meanTb_harmo'
    int_time = 6
    save_nc = True
    plot_ts_Tb_Tsys = False

    basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    
    #basename_lvl1 = "/scratch/GROSOM/Level1/"
    #basename_lvl2 = "/scratch/GROSOM/Level2/"

    if instrument_name=="GROMOS":
        import gromos_classes as gc
        calibration = gc.IntegrationGROMOS(date, basename_lvl1, integration_strategy, int_time)
    elif instrument_name=="SOMORA":
        import somora_classes as sm
        calibration = sm.IntegrationSOMORA(date, basename_lvl1, integration_strategy, int_time)
    elif instrument_name=="mopi5":
        import mopi5_classes as mc
        basename_lvl1 = "/scratch/MOPI5/Level1/"
        basename_lvl2 = "/scratch/MOPI5/Level2/"
        basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
        basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        #calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time, ['AC240','USRP-A'])
        calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time)

    calibrated_data, cal_flags, meteo_data = calibration.read_level1a()
    assert calibration.instrument_name == instrument_name, 'Wrong instrument definition'

    calibrated_data = calibration.clean_level1a_byFlags()

    calibrated_data = calibration.find_bad_channels_stdTb(spectrometers = calibration.spectrometers, stdTb_threshold = 10, apply_on='cal')

    #calibrated_data = calibration.add_mean_Tb(spectrometers = calibration.spectrometers)
    calibrated_data = calibration.add_mean_Tb(spectrometers = calibration.spectrometers, around_center=True, around_center_value=20e6)
    
    if plot_ts_Tb_Tsys:
        fig = plt.figure()
        ax = fig.add_subplot(2,1,1)
        ax1 = fig.add_subplot(2,1,2)
        #ax2 = fig.add_subplot(3,1,3)
        for s in calibration.spectrometers:
            ax.plot(calibration.calibrated_data[s].time, calibration.calibrated_data[s].mean_Tb,'.',markersize=4,label=s)
            ax1.plot(calibration.calibrated_data[s].time, calibration.calibrated_data[s].TSys,'.',markersize=4,label=s)
            #ax2.plot(calibration.calibration_flags[s].time, calibration.calibration_flags[s].sum(dim='flags').to_array()[0,:],'.',markersize=4,label=s)

            #ax.legend()        
            ax1.legend()    
        ax.set_ylabel('$T_B$ [K]')
        ax1.set_ylabel('$T_{sys}$ [K]')
        #plt.suptitle('Mean T_B and T_sys')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        #fig.savefig(calibration.level1_folder+'Tb_Tsys_all_'+calibration.datestr+'.pdf')

    # WARNING, depending on the integration type, some variable becomes meaningless --> for instance stdTb !!
    #integrated_data = calibration.integrate(spectrometers = calibration.spectrometers, strategy=integration_strategy, Tb_chunks=[150])

# %%
        
    # Define the parameters for integration
    TOD = [3, 9, 15, 21]
    interval = [3, 3, 3, 3]
    #meanTb_chunks = [80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 190]
    meanTb_chunks = [100, 110, 120, 130, 140, 160, 180]
    classic = np.arange(1,24)

    integrated_data, integrated_meteo = calibration.integrate(
        spectrometers = calibration.spectrometers, 
        strategy=integration_strategy,
        harmonized=True,
        tod=TOD, 
        interval=interval,
        Tb_chunks=meanTb_chunks,
        spectro_basis='U5303',
        freq=int_time,
        )

    if integration_strategy == 'meanTb' or integration_strategy == 'meanTb_harmo':
        dimension=['chunks','channel_idx']
    else:
        dimension=['time','channel_idx']

    # Note that the stdTb used here is just: mean(stdTb)/sqrt(N)
    integrated_data = calibration.find_bad_channels_stdTb(
        spectrometers = calibration.spectrometers, 
        stdTb_threshold = 15,
        apply_on='int',
        dimension=dimension
        )

    integrated_data = calibration.correct_troposphere(
        calibration.spectrometers, 
        dimension[0], 
        method='Ingold_v1', 
        basis_spectro='AC240')

    integrated_data = calibration.add_noise_level(calibration.spectrometers, max_diff_level=10)

    if instrument_name == 'mopi5':
        calibration.compare_spectra_mopi5(
            dim=dimension[0], 
            idx=np.arange(0,len(meanTb_chunks)+1), 
            #idx=[0,1,2,3],
            save_plot = True, 
            #identifier=TOD,
            identifier=meanTb_chunks+[300],
            with_corr = False
        )

        integrated_data = calibration.add_interpolated_spectra(
            spectrometers=['AC240', 'USRP-A', 'U5303'],
            use_basis='U5303',
            dim=dimension[0],              
            plot_diff=True)

        calibration.plot_time_min_comp()
    else:
        calibration.compare_Tb_chunks(dim=dimension[0], idx=[0,1,2,3], Tb_corr = True)

    if save_nc:
        calibration.save_dataset_level1b(
            spectrometers = calibration.spectrometers, 
            datasets=[integrated_data, integrated_meteo], 
            groups=['spectrometer1','meteo'], 
            extra='')


# %%
