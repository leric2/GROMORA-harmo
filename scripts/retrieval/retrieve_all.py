#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:37:52 2020

@author: eric

Main script for GROMOS and SOMORA long-term retrievals.

It only contains the operational parameters for GROMOS and SOMORA FFTS retrievals and is used to launch yearly retrievals from standard level 1.

For further retrievals option, the user should have a look at "retrieve.py".

Example:
    On birg, it works using the GROMORA_retrievals conda environment which needs to be activated with:
        $ conda activate GROMORA_retrievals

    Then, the user needs to go within the retrievals folder and the script can be launched with:
        $ python retrieve_all.py
"""
import sys

# Adding the paths to pyretrievals and retrieval folder
sys.path.insert(0, '/home/jobuser/scripts/instruments/gromos/GROMORA-harmo/scripts/retrieval/')
sys.path.insert(0, '/home/jobuser/scripts/instruments/gromos/GROMORA-harmo/scripts/pyretrievals/')

import datetime
import os
import time
from abc import ABC

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
import gc as garbage

from utils_GROSOM import save_single_pdf

# For ARTS, we need to specify some paths
load_dotenv('/opt/arts/.env.stockhorn-arts24')

ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

# Some important global paths that might vary based on the computer
# It always assumes that the data are separated in different folders for each years 
# within these basefolder.
GROMOS_L1_BASEFOLDER = '/storage/tub/instruments/gromos/level1/GROMORA/v2/'
GROMOS_L2_BASEFOLDER = '/storage/tub/instruments/gromos/level2/GROMORA/v2/'

SOMORA_L1_BASEFOLDER = '/storage/tub/instruments/somora/level1/v2/'
SOMORA_L2_BASEFOLDER = '/storage/tub/instruments/somora/level2/oper/'

def retrieve_day(date, instrument_name, integration_strategy='classic', int_time = 1):
    '''
    Function performing daily retrieval of GROMOS or SOMORA ozone profiles.

    '''
    # Dictionnary containing all EXTERNAL retrieval parameters 
    retrieval_param = dict()

    # Saving the ARTS paths into retrieval_param
    retrieval_param["ARTS_DATA_PATH"] = ARTS_DATA_PATH
    retrieval_param["ARTS_BUILD_PATH"] = ARTS_BUILD_PATH
    retrieval_param["ARTS_INCLUDE_PATH"] = ARTS_INCLUDE_PATH

    # Implementation of the instrument class depending on instrument name:
    if instrument_name=="GROMOS":
        import gromos_classes as gc
        basename_lvl1 = os.path.join(GROMOS_L1_BASEFOLDER,str(date.year))
        basename_lvl2 = os.path.join(GROMOS_L2_BASEFOLDER,str(date.year))
        instrument = gc.GROMOS_LvL2(
            date,
            basename_lvl1,
            basename_lvl2,
            integration_strategy,
            int_time,
            extra_base=''
            )
        retrieval_param['increased_var_factor'] = 1
    elif instrument_name=="SOMORA":
        basename_lvl1 = os.path.join(SOMORA_L1_BASEFOLDER,str(date.year))
        basename_lvl2 = os.path.join(SOMORA_L2_BASEFOLDER,str(date.year))
        import somora_classes as sm
        instrument = sm.SOMORA_LvL2(
            date=date,
            basename_lvl1=basename_lvl1,
            basename_lvl2=basename_lvl2,
            integration_strategy=integration_strategy,
            integration_time=int_time,
            extra_base=''
        )
        retrieval_param['increased_var_factor'] = 1
    
    # Reading of integrated (level 1) data:
    if integration_strategy == 'classic':
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
            no_flag=False, meta_data=True, extra_base='')
    else:
        raise NotImplementedError(
            'TODO, implement reading level1b in non classical cases !')

    # Some main parameters for the retrievals to perform
    # Type:
    retrieval_param["retrieval_type"] = 2
    # Retrieval quantities:
    retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'
    # Verbosity:
    retrieval_param['verbose'] = 1
    # Stopping after FM and plotting it:
    retrieval_param['FM_only'] = False
    retrieval_param['show_FM'] = False

    # The date:
    retrieval_param['date'] = date

    # Some plotting options (all off for operational routines)
    retrieval_param["plot_meteo_ds"] = False
    retrieval_param["show_f_grid"] = False
    retrieval_param['plot_opacities'] = False
    retrieval_param['plot_o3_apriori_covariance'] = False

    # Function to define the default retrieval_param dictionary.
    retrieval_param = instrument.define_retrieval_param(retrieval_param)

    # Quick test on instrument name
    assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'

    # Select the spectrometer (only AC240 supported at the moment)
    spectro = 'AC240'
    spectro_dataset = integrated_dataset[spectro]

    # Retrievel only the integration cycles with sufficient number of calibrated spectra.
    cycles=np.where(flags[spectro].calibration_flags.data[:,0]==1)[0]
    #cycles = [0,1,2]
    if len(cycles) ==0:
        return 0   

    figure_list = []
    
    counter = 0
    save_str='.pdf'
    for c in cycles:
        print('######################################################################################')
        retrieval_param["integration_cycle"] = c
        print('retrieving cycle : ',c)

        if ~np.isnan(integrated_meteo[spectro].air_pressure[c].data):
            retrieval_param['p_surface'] = integrated_meteo[spectro].air_pressure[c].data
        else:
            retrieval_param['p_surface'] = instrument.standard_air_pressure
        
        if ~np.isnan(integrated_meteo[spectro].air_temperature[c].data):
            retrieval_param['T_surface'] = integrated_meteo[spectro].air_temperature[c].data
        else:
            retrieval_param['T_surface'] = instrument.standard_air_temperature
            
        try:
            if retrieval_param["retrieval_type"] == 1:
                retrieval_param["surface_altitude"] = 1e3
                retrieval_param["observation_altitude"] =  15e3
                ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(
                    ac, 
                    spectro_dataset, 
                    retrieval_param, 
                    title = 'retrieval_trop_corr',
                    figure_list = figure_list
                    )
                level2_cycle = ac.get_level2_xarray()
                save_str = str(retrieval_param["integration_cycle"])+'_Perrin_corr.pdf'
                #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
                #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
            elif retrieval_param["retrieval_type"] == 2:
                retrieval_param["surface_altitude"] = 1000
                retrieval_param["observation_altitude"] =  1000   
                ac, retrieval_param, sensor_out = instrument.retrieve_cycle(spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
                if ac.oem_converged:
                    #figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='ozone retrieval cycle' + str(c),figure_list=figure_list)
                    level2_cycle = ac.get_level2_xarray()
                else:
                    level2_cycle=xr.Dataset()
                save_str = str(retrieval_param["integration_cycle"])+'_all_retrieved_Perrin.pdf'
            else:
                save_str = '.pdf'
                level2_cycle=[]
            
            if (counter == 0) & (len(level2_cycle) > 0):
                counter = counter + 1
                level2 = level2_cycle
            elif (counter>0) & (len(level2_cycle) > 0):
                counter = counter + 1
                level2 = xr.concat([level2, level2_cycle], dim='time')
        except Exception as e:
            print(e)
            # if counter==0:
            #     level2 = xr.Dataset()
            # else:
            #     level2_cycle=xr.Dataset()
            print('problem retrieving cycle : ',c)

    if counter > 0:
        #save_single_pdf(instrument.filename_level2[spectro]+'_'+save_str, figure_list)
        level2 = instrument.write_level2_gromora(level2, retrieval_param, full_name = instrument.filename_level2[spectro]+'_v2.nc')
        level2.close()
        level2_cycle.close()
        del level2, level2_cycle
        garbage.collect()

if __name__ == "__main__":
    # Selection of the integration stategy used for the level 1. 
    # Currently only 'classic' supported for GROMOS and SOMORA.
    integration_strategy = 'classic'

    # A selection of days that crash the retrievals for unknown reasons. Not needed in latest version
    void_date_problem = []

    # Date range on which to perform the retrievals
    #dates = pd.date_range(start='2022-05-25', end='2022-05-25')#.append(pd.date_range(start='2010-01-01', end='2010-01-03')).append(pd.date_range(start='2015-01-01', end='2015-01-04'))#.append(pd.date_range(start='2012-11-26', end='2012-12-31'))#).append(pd.date_range(start='2016-12-31', end='2017-01-01'))
    date = pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=4))
    print('######################################################################################')
    print('######################################################################################')
    print('######################################################################################')

    
    #    if d in void_date_problem :
    #        print('abort core problem with this day : ',d ,' --> skipping')
    #    else:
            # try:
            #     retrieve_day(d, instrument_name='SOMORA', integration_strategy=integration_strategy)
            # except:
            #    pass
            # print('######################################################################################')
    print('######################################################################################')
   # try:
    retrieve_day(date, instrument_name='GROMOS',  integration_strategy=integration_strategy)
   # except:
   #     pass            

