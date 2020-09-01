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
    instrument_name = "mopi5"
    #date = datetime.date(2019,2,1)
    date = pd.date_range(start='2019-02-01', end='2019-02-01')

    integration_strategy = 'simple'
    int_time = 1

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
        #basename_lvl1 = "/scratch/MOPI5/Level1/"
        #basename_lvl2 = "/scratch/MOPI5/Level2/"
        basename_lvl1 = "/home/eric/Documents/PhD/DATA/"
        basename_lvl2 = "/home/eric/Documents/PhD/DATA/"
        calibration = mc.IntegrationMOPI5(date, basename_lvl1, integration_strategy, int_time)

    calibrated_data, cal_flags, cal_meteo = calibration.read_level1a()
    assert calibration.instrument_name == instrument_name, 'Wrong instrument definition'

    calibration.integrate(calibration.spectrometers)