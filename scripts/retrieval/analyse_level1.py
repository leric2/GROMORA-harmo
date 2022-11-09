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
import datetime
import os
from abc import ABC

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import scipy.io
import xarray as xr
from dotenv import load_dotenv


plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]})
load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')

instrument_name = "GROMOS"

# date = pd.date_range(start='2019-01-03', end='2019-01-05')
# meanTb_chunks = [95, 100, 110, 120, 130, 140, 180]
# lowerBound = [0, 95, 100, 110, 120, 130, 140, 180]

# date = pd.date_range(start='2019-01-30', end='2019-06-18')

#date = pd.date_range(start='2018-01-10', end='2018-01-11')
#date = pd.date_range(start='2017-09-01', end='2018-01-05')
date = datetime.date(2019,10,3)
#date = [datetime.date(2019,3,11), datetime.date(2019,4,3)]

int_time = 1

if instrument_name == "GROMOS":
    import gromos_classes as gc
    basename_lvl1 = "/storage/tub/instruments/gromos/level1/GROMORA/"+str(date.year)
    #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v1/"+str(date.year)
    instrument = gc.GROMOS_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy='classic',
        integration_time=int_time
    )
elif instrument_name == "SOMORA":
    import somora_classes as sm
    basename_lvl1 = "/storage/tub/instruments/somora/level1/v1/"+str(date.year)
    #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
    basename_lvl2 = "/storage/tub/instruments/somora/level2/v1/"+str(date.year)
    instrument = sm.SOMORA_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy='classic',
        integration_time=int_time
    )

def read_gromos_old_level1(filename, cycle=0):
    basename = '/home/esauvageat/Documents/OG_GROMOS/meta3/level1/'
    full_name = os.path.join(basename, filename)
    old_level1 = scipy.io.loadmat(full_name)
    old_level1 = old_level1['level1_data']

    Tb = old_level1[0,cycle][3]
    freq = old_level1[0,cycle][4]

    return [freq, Tb]

def read_somora_old_level1(filename_somora, cycle=0):
    basename = '/home/esauvageat/Documents/SOMORA/Level1/'
    full_name = os.path.join(basename, filename_somora)
    old_level1 = scipy.io.loadmat(full_name)
    old_level1 = old_level1['level1_data']

    Tb = old_level1[0,cycle][3]

    center=1.4217504e+11
    somora_frequencies=center-(8191*(1e+9/16384)):(1e+9/16384):center+(8192*(1e+9/16384))

    return [freq, Tb]

def compare_old_new_tb(old_freq, old_Tb, new, cycle):
    s = 'AC240'
    fig, axs = plt.subplots(1,1,sharex=True)
    
    colors =  ['#d95f02','#1b9e77']
    mask = new[s].good_channels[cycle].data
    mask[mask==0]=np.nan
    axs.plot(new[s].frequencies[cycle].data/1e9,new[s].Tb[cycle]*mask, lw=0.5, color='r', label='new routine')
    axs.plot(old_freq/1e9,old_Tb, lw=0.5, color='b', label='old routine')

    #axs.set_xlim(110.25, 111.4)
    axs.set_ylim(80,120)
    #ax].set_ylim(np.median(ds_dict[s].Tb[calibration_cycle].data)-10,np.median(ds_dict[s].Tb[calibration_cycle].data)+15)
    axs.set_xlabel("frequency [GHz]")
    axs.set_ylabel(r"$T_B$ [K]")
        
    axs.set_title('')
    axs.grid()
    #cbar.ax.set_yticklabels(np.arange(0,100,1))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# Reading data
level1b_dataset = instrument.read_level1b(no_flag=True)
F0 = instrument.observation_frequency

filename_somora = 'm20190320.bin'


filename = 'level1_dt_60min_2019_10_03.mat'
old_freq, old_Tb = read_gromos_old_level1(filename, cycle=12)

compare_old_new_tb(old_freq, old_Tb, level1b_dataset[0], 12)

# Plotting selected spectra
instrument.plot_level1b_TB(title=instrument_name+' ',save=True, outfolder='/scratch/GROSOM/Level1/',save_name=instrument_name+'_integrated_spectra_',idx=[12])

