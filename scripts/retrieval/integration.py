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

# %%

def integrate(date, integration_strategy):
    instrument_name = "mopi5"
    #date = datetime.date(2019,2,21)
    #date = pd.date_range(start='2019-01-03', end='2019-01-05')
    meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]


    #date = pd.date_range(start='2019-01-30', end='2019-06-18')

    #date = pd.date_range(start='2019-01-30', end='2019-02-22')
    #meanTb_chunks = [80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 170, 190]

    #date = pd.date_range(start='2019-05-01', end='2019-05-04')
    # No U5303

    #date = pd.date_range(start='2019-04-25', end='2019-04-27')
    #meanTb_chunks = [105, 110, 115, 120, 130, 160, 180, 200]


    #date = pd.date_range(start='2019-06-11', end='2019-06-15')
    #meanTb_chunks = [110, 120, 130, 140, 150, 160, 170, 180, 200, 220]


    #date = pd.date_range(start='2019-03-12', end='2019-03-12')
    #date = pd.date_range(start='2019-06-13', end='2019-06-13')
    # options are: 'TOD', 'TOD_harmo', 'classic' 'meanTb_harmo', or 'meanTb'
    #integration_strategy = 'meanTb_harmo'
    int_time = 1
    save_nc = True
    plot_ts_Tb_Tsys = False
    df_bins=200e3

    basename_lvl1 = "/home/eric/Documents/PhD/DATA/Level1a/"
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
        basename_lvl1 = "/home/eric/Documents/PhD/DATA/MOPI5/"
        basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        #calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time, ['AC240','USRP-A'])
        calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time)

    # Define the parameters for integration
    #TOD = [3, 9, 15, 21]
    #interval = [3, 3, 3, 3]
    #TOD = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]
    #interval = np.ones(len(TOD))
    TOD = np.arange(24)
    interval = 0.5*np.ones(len(TOD))
    
    #meanTb_chunks=[95, 100, 105]
    
    classic = np.arange(1,24)

    if integration_strategy == 'meanTb' or integration_strategy == 'meanTb_harmo':
        dimension=['chunks','channel_idx']
    else:
        dimension=['time','channel_idx']

# %%
    calibrated_data, cal_flags, meteo_data = calibration.read_level1a()
    assert calibration.instrument_name == instrument_name, 'Wrong instrument definition'

    calibrated_data = calibration.clean_level1a_byFlags_all()

    calibrated_data = calibration.find_bad_channels_stdTb(spectrometers = calibration.spectrometers, stdTb_threshold = 10, apply_on='cal')

    #calibrated_data = calibration.add_mean_Tb(spectrometers = calibration.spectrometers)
    calibrated_data = calibration.add_mean_Tb(spectrometers = calibration.spectrometers, around_center=True, around_center_value=20e6)
    
    if plot_ts_Tb_Tsys:
        calibration.plot_time_series_all_mopi5(special=False)
        #fig.savefig(calibration.level1_folder+'Tb_Tsys_all_'+calibration.datestr+'.pdf')

    # WARNING, depending on the integration type, some variable becomes meaningless --> for instance stdTb !!
    #integrated_data = calibration.integrate(spectrometers = calibration.spectrometers, strategy=integration_strategy, Tb_chunks=[150])

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

    # Note that the stdTb used here is just: mean(stdTb)/sqrt(N)
    integrated_data = calibration.find_bad_channels_stdTb(
        spectrometers = calibration.spectrometers, 
        stdTb_threshold = 15,
        apply_on='int',
        dimension=dimension
        )

    # skip_ch_U5303 = [int(integrated_data['U5303'].channel_idx.where(integrated_data['U5303'].frequencies[3]>calibration.observation_frequency-400e6, drop=True)[0].data),
    #                 len(calibration.integrated_data['U5303'].channel_idx)-int(integrated_data['U5303'].channel_idx.where(integrated_data['U5303'].frequencies[3]>calibration.observation_frequency+400e6, drop=True)[3].data)]
    # skip_ch_AC240 = [int(integrated_data['AC240'].channel_idx.where(integrated_data['AC240'].frequencies[3]>calibration.observation_frequency-400e6, drop=True)[0].data),
    #                 len(calibration.integrated_data['AC240'].channel_idx)-int(integrated_data['AC240'].channel_idx.where(integrated_data['AC240'].frequencies[3]>calibration.observation_frequency+400e6, drop=True)[3].data)]
    
    skip_ch_U5303 = [1541, 6651]
    skip_ch_AC240 = [2559, 815]
    integrated_data = calibration.correct_troposphere(
        ['U5303', 'USRP-A'], 
        dimension[0], 
        method='Ingold_v1', 
        basis_spectro='U5303',
        skip_ch = skip_ch_U5303,
        num_of_ch = 500
        )
        
    integrated_data = calibration.correct_troposphere(
        ['AC240'], 
        dimension[0], 
        method='Ingold_v1', 
        basis_spectro='AC240',
        skip_ch = skip_ch_AC240,
        num_of_ch = 800
        )

    integrated_data = calibration.add_noise_level(calibration.spectrometers, max_diff_level=10)
    
    if instrument_name == 'mopi5':
        
        integrated_data = calibration.add_binned_spectra(
            spectrometers=['AC240', 'USRP-A', 'U5303'],
            bin_size=[3, 15, 2],
            dim=dimension[0],
            df=None,              
            plot_diff=False)
        
        integrated_data = calibration.add_binned_spectra(
            spectrometers=['AC240', 'USRP-A', 'U5303'],
            bin_size=[3, 15, 2],
            dim=dimension[0],
            df=None,              
            plot_diff=False,
            corrected=True
            )

        integrated_data = calibration.add_interpolated_spectra(
            spectrometers=['AC240', 'USRP-A', 'U5303'],
            use_basis='U5303',
            dim=dimension[0],              
            plot_diff=False,
            from_binned=True
            )

        integrated_data = calibration.add_interpolated_spectra(
            spectrometers=['AC240', 'USRP-A', 'U5303'],
            use_basis='U5303',
            dim=dimension[0],              
            plot_diff=False,
            from_binned=True,
            corrected=True
            )
        
        # integrated_data = calibration.correct_troposphere(
        #     ['U5303', 'USRP-A'], 
        #     dimension[0], 
        #     method='Ingold_v1', 
        #     basis_spectro='U5303',
        #     skip_ch = [1000,5300],
        #     num_of_ch = 500,
        #     interp=True
        #     )

        # integrated_data = calibration.correct_troposphere(
        #     ['AC240'], 
        #     dimension[0], 
        #     method='Ingold_v1', 
        #     basis_spectro='AC240',
        #     skip_ch = [1000,1000],
        #     num_of_ch = 500,
        #     interp=True
        #     )
        #param_slope_broadband = [111.2e9, 25e6, 110.5e9, 25e6]
        param_slope_broadband = [111.036e9, 25e6, 110.336e9, 25e6]
        param_slope = {'AC240' : param_slope_broadband, 'USRP-A': [110.84e9, 5e6, 110.72e9, 5e6], 'U5303':param_slope_broadband}
        integrated_data = calibration.add_bias_characterization_new(
            spectrometers=calibration.spectrometers,
            dim=dimension[0],
            param_slope = param_slope, 
            around_center_value=1e6
            )
        # integrated_data = calibration.add_bias_characterization(
        #     spectrometers=['AC240', 'USRP-A'],
        #     use_basis='U5303',
        #     dim=dimension[0],
        #     param_slope = param_slope, 
        #     around_center_value=1e6
        #     )

        integrated_data = calibration.add_bias_characterization_corrected(
            spectrometers=['AC240', 'USRP-A'],
            use_basis='U5303',
            dim=dimension[0],
            param_slope = param_slope, 
            around_center_value=1e6
            )

    if save_nc:
        calibration.save_dataset_level1b(
            spectrometers = calibration.spectrometers, 
            datasets=[integrated_data, integrated_meteo], 
            groups=['spectrometer1','meteo'], 
            extra='')

# %%

def plot_integrated(date, integration_strategy):
    #date1b = pd.to_datetime(date[-1])
    date1b = pd.to_datetime(date[-1])
    if instrument_name=="GROMOS":
        import gromos_classes as gc
        integration = gc.GROMOS_LvL2(date, basename_lvl1, integration_strategy, int_time)
    elif instrument_name=="SOMORA":
        import somora_classes as sm
        integration = sm.SOMORA_LvL2(date, basename_lvl1, integration_strategy, int_time)
    elif instrument_name=="mopi5":
        import mopi5_classes as mc
        basename_lvl1 = "/scratch/MOPI5/Level1/"
        basename_lvl2 = "/scratch/MOPI5/Level2/"
        #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
        #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        
        integration = mc.MOPI5_LvL2(date1b, basename_lvl1, basename_lvl2, integration_strategy, integration_time=int_time)

    # Plotting part
    integrated_data, integrated_flags, integrated_meteo = integration.read_level1b(no_flag=True, meta_data=False)

    if instrument_name == 'mopi5':
        integration.compare_spectra_mopi5(
            dim=dimension[0], 
            idx=np.arange(0,len(meanTb_chunks)+1), 
            #idx=np.arange(len(TOD)),
            save_plot = True, 
            #identifier=TOD,
            identifier=meanTb_chunks+[300],
            with_corr = True
        )
        integration.plot_time_min_comp()
    
        integration.compare_spectra_binned_interp_mopi5(
                dim=dimension[0], 
                idx=np.arange(0,len(meanTb_chunks)+1), 
                #idx=np.arange(len(TOD)),
                spectrometers=integration.spectrometers,
                save_plot = True, 
                #identifier=TOD,
                identifier=meanTb_chunks+[300],
                use_basis='U5303',
                corrected=True
        )
        # integration.compare_binned_spectra_mopi5(
        #         dim=dimension[0], 
        #         idx=np.arange(0,len(meanTb_chunks)+1), 
        #         #idx=[0,1,2,3],
        #         save_plot = True, 
        #         #identifier=TOD,
        #         identifier=meanTb_chunks+[300],
        #         use_basis='U5303',
        #         df=df_bins
        #     )
        
        # integration.compare_interpolated_spectra_mopi5(
        #     dim=dimension[0], 
        #     idx=np.arange(0,len(meanTb_chunks)+1), 
        #     #idx=[0,1,2,3],
        #     spectrometers=integration.spectrometers,
        #     use_basis='U5303',
        #     save_plot = True, 
        #     #identifier=TOD,
        #     identifier=meanTb_chunks+[300],
        # )

    else:
        integration.compare_Tb_chunks(dim=dimension[0], idx=[0,1,2,3], Tb_corr = True)

# %%
if __name__ == "__main__":
    dateR = pd.date_range(start='2019-01-03', end='2019-01-05')
    #dateR = pd.date_range(start='2019-04-25', end='2019-04-27')
    #dateR = pd.date_range(start='2019-01-30', end='2019-02-22')
    integrate(dateR, 'meanTb_harmo')
    
    # options are: 'TOD', 'TOD_harmo', 'classic' 'meanTb_harmo', or 'meanTb'
    # integration_strategy = 'TOD_harmo'

    # for date in dateR:
    #      try:
    #          integrate(date, integration_strategy)
    #      except:
    #          print('not working for day : ', date)
        #plot_integrated(date, integration_strategy)

# %%

