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
import sys

sys.path.insert(0, '/home/es19m597/Documents/GROMORA/GROMORA-harmo/scripts/retrieval/')
sys.path.insert(0, '/home/es19m597/Documents/GROMORA/GROMORA-harmo/scripts/pyretrievals/')


import datetime
import os
import time
from abc import ABC
import sys

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
import gc

from gromora_utils import save_single_pdf

sys.path.insert(0, '/home/esauvageat/Documents/GROMORA/Analysis/GROMORA-harmo/scripts/retrieval/')
sys.path.insert(0, '/home/esauvageat/Documents/GROMORA/Analysis/GROMORA-harmo/scripts/pyretrievals/')


# For ARTS, we need to specify some paths
load_dotenv('/opt/anaconda/.env.birg-arts24_pyarts')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

def sensitivity_analysis(instrument_name, date, param, cycles):
    for test in param:
        int_time = 1
        integration_strategy = 'classic'
        recheck_channels = False

        basename_lvl2 = "/home/es19m597/Documents/GROMORA/Data/"

        # Dictionnary containing all EXTERNAL retrieval parameters 
        retrieval_param = dict()
        retrieval_param["ARTS_DATA_PATH"] = ARTS_DATA_PATH
        retrieval_param["ARTS_BUILD_PATH"] = ARTS_BUILD_PATH
        retrieval_param["ARTS_INCLUDE_PATH"] = ARTS_INCLUDE_PATH

        if instrument_name == "GROMOS":
            import gromos_classes as gc
            basename_lvl1 = os.path.join(
                '/storage/tub/instruments/gromos/level1/GROMORA/v2/', str(date.year))
            instrument = gc.GROMOS_LvL2(
                date,
                basename_lvl1,
                basename_lvl2,
                integration_strategy,
                int_time,
                extra_base=''
                )
            retrieval_param['increased_var_factor'] = 1  # 15
        elif instrument_name == "SOMORA":
            basename_lvl1 = os.path.join(
                '/storage/tub/instruments/somora/level1/v2/', str(date.year))
            import somora_classes as sm
            instrument = sm.SOMORA_LvL2(
                date=date,
                basename_lvl1=basename_lvl1,
                basename_lvl2=basename_lvl2,
                integration_strategy=integration_strategy,
                integration_time=int_time,
                extra_base=''
            )
            retrieval_param['increased_var_factor'] = 1  # 1.1 #15

        if integration_strategy == 'classic':
            integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
                no_flag=False, meta_data=True, extra_base='')
        else:
            raise NotImplementedError(
                'TODO, implement reading level1b in non classical cases !')
        retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'
        retrieval_param['verbose'] = 2
        retrieval_param["retrieval_type"] = 8
        retrieval_param['FM_only'] = False
        retrieval_param['show_FM'] = True
        retrieval_param['sensor'] = 'FFT_SB'
        retrieval_param['SB_bias'] = 0
        retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'

        retrieval_param["plot_meteo_ds"] = False

        retrieval_param["show_f_grid"] = False
        retrieval_param['plot_opacities'] = False
        # 3. test retrieving the FM
        retrieval_param['verbose'] = 1

        retrieval_param['date'] = date

        retrieval_param['plot_o3_apriori_covariance'] = True


        retrieval_param = instrument.define_retrieval_param(retrieval_param)

        assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'

        if recheck_channels:
            integrated_data = instrument.find_bad_channels_stdTb(
                spectrometers=instrument.spectrometers,
                stdTb_threshold=12,
                apply_on='int',
                dimension=['time', 'channel_idx']
            )

        retrieval_param['plot_o3_apriori_covariance'] = False

        retrieval_param = instrument.define_retrieval_param(retrieval_param)

        # for i,s in enumerate(instrument.spectrometer):
        spectro = 'AC240'
        spectro_dataset = instrument.integrated_data[spectro]

        figure_list = []

        counter = 0
        for c in cycles:
            counter = counter + 1
            retrieval_param["integration_cycle"] = c
            print('retrieving cycle : ', c)

            retrieval_param['p_surface'] = integrated_meteo[spectro].air_pressure[c].data
            retrieval_param['T_surface'] = integrated_meteo[spectro].air_temperature[c].data

            if retrieval_param["retrieval_type"] == 8:
                retrieval_param["surface_altitude"] = 1000
                retrieval_param["observation_altitude"] = 1000
                retrieval_param["test_type"] = test
                print('Sensitivity Analysis in :', test)
                if test=='og':
                    save_str = 'sensitivity_test_og'
                elif test == 'angle':
                    save_str = 'sensitivity_test_angle'
                    retrieval_param['pointing_angle_corr'] = 1
                elif test=='neg_angle':
                    save_str = 'sensitivity_test_angle_neg'
                    retrieval_param['pointing_angle_corr'] = -1
                elif test=='Tprofile':
                    save_str = 'sensitivity_test_Tprofile'
                    retrieval_param['Tprofile_bias'] = 5
                elif test=='continuum':
                    # 'H2O-PWR98' #H2O-ContMPM93 'H2O-PWR98, H2O' #"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
                    retrieval_param['water_vapor_model'] = 'H2O-MPM93'
                    retrieval_param['selected_species'] = ['O3', retrieval_param['water_vapor_model'],
                        retrieval_param['o2_model'], retrieval_param['n2_model']]
                    save_str='sensitivity_test_continuum'
                elif test=='noise':
                    retrieval_param['increased_var_factor'] = 1.1
                    save_str='sensitivity_test_noise'
                elif test=='SB':
                    retrieval_param['SB_bias'] = +0.05e-3 # 0.05mm uncertainties, total is 11mm for SOMORA and 20mm for GROMOS
                    save_str = 'sensitivity_test_SB'
                elif test=='Tcold' or test=='tWindow':
                    extra = '_'+test
                    save_str=''
                    #save_str = test

                    if instrument_name == "GROMOS":
                        instrument = gc.GROMOS_LvL2(
                            date,
                            basename_lvl1,
                            basename_lvl2,
                            integration_strategy,
                            int_time,
                            extra_base=extra
                            )
                    elif instrument_name == "SOMORA":
                        instrument = sm.SOMORA_LvL2(
                            date=date,
                            basename_lvl1=basename_lvl1,
                            basename_lvl2=basename_lvl2,
                            integration_strategy=integration_strategy,
                            integration_time=int_time,
                            extra_base=extra
                        )
                    integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
                        no_flag=False, meta_data=True, extra_base=extra)
                    spectro_dataset = instrument.integrated_data[spectro]
                elif test=='spectroscopy':
                    save_str = 'sensitivity_test_spectroscopy_new'
                    #line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
                    #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"
                    #line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all_species.par'
                    retrieval_param['FM_only'] = True
                    line_file = '/home/es19m597/Documents/GROMORA/InputsRetrievals/Perrin_modified_linestrength.xml'
                    line_file = '/home/es19m597/Documents/GROMORA/InputsRetrievals/Perrin_modified.xml'
                    retrieval_param['line_file'] = line_file
                else: 
                    save_str='other'

                ac, retrieval_param, sensor_out = instrument.retrieve_cycle(spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)

                figure_list = instrument.plot_level2(
                    ac, spectro_dataset, retrieval_param, title='retrieval_o3', figure_list=figure_list)
                level2_cycle = ac.get_level2_xarray()
                save_single_pdf(
                    instrument.filename_level2[spectro]+'_'+save_str+'.pdf', figure_list)
                level2_cycle.to_netcdf(
                    path=instrument.filename_level2[spectro]+'_'+save_str+'.nc')
            else:
                level2_cycle = xr.Dataset()

if __name__=='__main__':
    instrument_name = ['GROMOS'] #,'SOMORA'
    #date = datetime.date(2018, 6, 9)
    date = datetime.date(2018, 2, 26)
    #cycle = np.arange(7, 8)
    cycle = np.arange(9, 10)

    tests = ['spectroscopy'] # 'og', 'continuum','angle','spectroscopy','Tprofile','SB','Tcold','tWindow']
    for radiometer in instrument_name:
        try:
            sensitivity_analysis(radiometer, date, param=tests, cycles=cycle)
        except:
            pass