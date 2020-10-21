#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Retrieval script

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
    instrument_name = "mopi5"
    date = datetime.date(2019,1,5)
    int_time = 1
    integration_strategy = 'classic'
    recheck_channels = False
    
    basename_lvl1 = "/scratch/GROSOM/Level1/"
    basename_lvl2 = "/scratch/GROSOM/Level2/"

    if instrument_name=="GROMOS":
        import gromos_classes as gc
        instrument = gc.GROMOS_LvL2(
            date, 
            basename_lvl1, 
            basename_lvl2, 
            integration_strategy, 
            int_time)
    elif instrument_name=="SOMORA":
        import somora_classes as sm
        instrument = sm.SOMORA_LvL2(
            date=date,
            basename_lvl1=basename_lvl1,
            basename_lvl2=basename_lvl2,
            integration_strategy=integration_strategy,
            integration_time=int_time
        )
    elif instrument_name=="mopi5":
        import mopi5_classes as mc
        basename_lvl1 = "/scratch/MOPI5/Level1/"
        basename_lvl2 = "/scratch/MOPI5/Level2/"
        #basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
        #basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        instrument = mc.MOPI5_LvL2(
            date=date,
            basename_lvl1=basename_lvl1,
            basename_lvl2=basename_lvl2,
            integration_strategy=integration_strategy,
            integration_time=int_time
        )
    
    # Dictionnary containing all EXTERNAL retrieval parameters 
    retrieval_param = dict()

    if integration_strategy == 'classic':
        integrated_data, flags, integrated_meteo = instrument.read_level1b()
    else:
        raise NotImplementedError('TODO, implement reading level1b in non classical cases !')

    assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'

    if recheck_channels:
        integrated_data = instrument.find_bad_channels_stdTb(
            spectrometers = instrument.spectrometers, 
            stdTb_threshold = 12,
            apply_on='int',
            dimension=['time','channel_idx']
            )

    base_title = "GROMOS hourly integrated spectra : "
    #instrument.plot_level1b_TB_all(title = base_title, save=True, save_name='integrated_spectra_all_')
    base_title = "GROMOS hourly integrated spectrum : "
    #instrument.plot_level1b_TB_all(title = base_title, save=True, save_name='integrated_spectra_sel_',idx=[0])