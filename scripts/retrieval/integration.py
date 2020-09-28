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

if __name__ == "__main__":
    instrument_name = "GROMOS"
    #date = datetime.date(2019,2,1)
    date = pd.date_range(start='2019-01-30', end='2019-02-22')
    #date = pd.date_range(start='2019-05-01', end='2019-05-04')
    #date = pd.date_range(start='2019-04-25', end='2019-04-27')
    #date = pd.date_range(start='2019-06-11', end='2019-06-18')
    date = pd.date_range(start='2019-02-10', end='2019-02-10')
    # options are: 'TOD', 'classic' or 'meanTb'
    integration_strategy = 'classic'
    int_time = 1
    save_nc = False

    #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
    #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
    
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/"

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
        #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
        #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        #calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time, ['AC240','USRP-A'])
        calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time)

    calibrated_data, cal_flags, meteo_data = calibration.read_level1a()
    assert calibration.instrument_name == instrument_name, 'Wrong instrument definition'

    calibrated_data = calibration.clean_level1a_byFlags()

    calibrated_data = calibration.find_bad_channels_stdTb(spectrometers = calibration.spectrometers, stdTb_threshold = 10, apply_on='cal')

    #calibrated_data = calibration.add_mean_Tb(spectrometers = calibration.spectrometers)
    calibrated_data = calibration.add_mean_Tb(spectrometers = calibration.spectrometers, around_center=True, around_center_value=50e6)
    
    # WARNING, depending on the integration type, some variable becomes meaningless --> for instance stdTb !!
    #integrated_data = calibration.integrate(spectrometers = calibration.spectrometers, strategy=integration_strategy, Tb_chunks=[150])
    
    # Define the parameters for integration
    TOD = [4, 22]
    interval = [0.5, 0.5]
    meanTb_chunks = [100, 125, 150, 175]

    integrated_data, integrated_meteo = calibration.integrate(
        spectrometers = calibration.spectrometers, 
        strategy=integration_strategy, 
        tod=TOD, 
        interval=interval,
        Tb_chunks=meanTb_chunks,
        freq=int_time,
        )

    if integration_strategy == 'meanTb':
        dimension=['chunks','channel_idx']
    else:
        dimension=['time','channel_idx']

    # Note that the stdTb used here is just: mean(stdTb)/sqrt(N)
    integrated_data = calibration.find_bad_channels_stdTb(
        spectrometers = calibration.spectrometers, 
        stdTb_threshold = 8,
        apply_on='int',
        dimension=dimension
        )

    integrated_data = calibration.correct_troposphere(calibration.spectrometers, dimension[0])

    integrated_data = calibration.add_noise_level(calibration.spectrometers, max_diff_level=10)

    if instrument_name == 'mopi5':
        calibration.compare_spectra_mopi5(dim=dimension[0], idx=[0,1,2,3], save_plot = True, identifier=meanTb_chunks)
    else:
        calibration.compare_Tb_chunks(dim=dimension[0], idx=[0,1,2,3], Tb_corr = True)

    if save_nc:
        calibration.save_dataset_level1b(
            spectrometers = calibration.spectrometers, 
            datasets=[integrated_data, integrated_meteo], 
            groups=['spectrometer1','meteo'], 
            extra='')

    