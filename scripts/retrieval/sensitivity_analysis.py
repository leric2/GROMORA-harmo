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

from utils_GROSOM import save_single_pdf

# For ARTS, we need to specify some paths
load_dotenv('/opt/anaconda/.env.birg-arts24')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

def sensitivity_analysis(instrument_name, date, param, cycles):
    for test in param:
        int_time = 1
        integration_strategy = 'classic'
        recheck_channels = False

        basename_lvl2 = "/home/es19m597/Documents/GROMORA/Data/"

        line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
        #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"
        #line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all_species.par'
        #line_file = '/home/es19m597/Documents/GROMORA/InputsRetrievals/Hitran_f_modified.xml'

        # Dictionnary containing all EXTERNAL retrieval parameters
        retrieval_param = dict()

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
        elif instrument_name == "mopi5":
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

        if integration_strategy == 'classic':
            integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
                no_flag=False, meta_data=True, extra_base='')
        else:
            raise NotImplementedError(
                'TODO, implement reading level1b in non classical cases !')


        # type of retrieval to do:
        # 1. tropospheric corrected
        # 2. with h20
        # 3. test retrieving the FM
        retrieval_param["retrieval_type"] = 8
        retrieval_param['FM_only'] = False
        retrieval_param['show_FM'] = True
        retrieval_param['sensor'] = 'FFT_SB'
        retrieval_param['SB_bias'] = 0
        retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit'

        retrieval_param["obs_freq"] = instrument.observation_frequency
        retrieval_param['sideband_response'] = 'theory'
        retrieval_param["plot_meteo_ds"] = False

        retrieval_param["show_f_grid"] = True
        retrieval_param['plot_opacities'] = False
        retrieval_param["f_shift"] = 0  # 70e3#-405e3# -324.992e3#+1500e3
        retrieval_param["number_of_freq_points"] = 1201
        retrieval_param["irregularity_f_grid"] = 45
        retrieval_param["z_top_sim_grid"] = 112e3
        retrieval_param["z_bottom_sim_grid"] = 600
        retrieval_param["z_resolution_sim_grid"] = 2e3

        retrieval_param["retrieval_grid_type"] = 'altitude'
        retrieval_param["z_top_ret_grid"] = 95e3
        retrieval_param["z_bottom_ret_grid"] = 1e3
        retrieval_param["z_resolution_ret_grid"] = 2e3

        retrieval_param["retrieval_h2o_grid_type"] = 'pressure'
        # retrieval_param["z_top_ret_grid_h2o"] = 20e3
        # retrieval_param["z_bottom_ret_grid_h2o"] = 1e3
        # retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

        retrieval_param["h2o_pressure"] = [500e2]

        #retrieval_param['unit_var_y']  = 3**2
        retrieval_param['pointing_angle_corr'] = 0

        retrieval_param['apriori_ozone_climatology_GROMOS'] = '/storage/tub/instruments/gromos/InputsRetrievals/apriori_ECMWF_MLS/'
        retrieval_param['apriori_ozone_climatology_SOMORA'] = '/storage/tub/instruments/gromos/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'

        #retrieval_param['obs_freq'] = 1.4217504e11
        retrieval_param['spectroscopy_type'] = 'XML'
        retrieval_param['line_file'] = line_file
        retrieval_param['atm'] = 'ecmwf_cira86'  # fascod  ecmwf_cira86
        # max_diff, simple_stack_corr, simple
        retrieval_param['ptz_merge_method'] = 'max_diff'
        retrieval_param['ptz_merge_max_Tdiff'] = 5
        retrieval_param['h2o_apriori'] = 'ecmwf'  # 'ecmwf' # 'fascod_extended'
        # /tub/instruments/gromos/ECMWF_Bern'
        retrieval_param['ecmwf_store_location'] = '/storage/tub/atmosphere/ecmwf/locations/'+str(date.year)
        #retrieval_param['ecmwf_store_location'] ='/home/eric/Documents/PhD/ECMWF'
        retrieval_param['extra_time_ecmwf'] = 3.5

        retrieval_param['o3_apriori'] = 'waccm_monthly'#'waccm_monthly'
        # 'waccm_yearly_scaled'low_alt_ratio
        retrieval_param['plot_o3_apriori_covariance'] = True
        retrieval_param['o3_apriori_covariance'] = 'low_alt_ratio_optimized' #low_alt_ratio_optimized
        retrieval_param['waccm_file'] = '/storage/tub/instruments/gromos/InputsRetrievals/waccm_o3_climatology.nc'
        retrieval_param["apriori_O3_cov"] = 1e-6  # 1e-6
        retrieval_param["apriori_H2O_stdDev"] = 1  # 6e-4 12e-4 0.5 16e-4

        #retrieval_param["apriori_o2_stdDev"] = 1e-8  # 6e-4
        #retrieval_param["apriori_n2_stdDev"] = 1e-8

        retrieval_param['covmat_polyfit_0'] = 0.1
        retrieval_param['covmat_polyfit_1'] = 0.5
        retrieval_param['covmat_polyfit_2'] = 0.5

        # 'H2O-PWR98' #H2O-ContMPM93 'H2O-PWR98, H2O' #"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
        retrieval_param['water_vapor_model'] = 'H2O-PWR98'
        retrieval_param['o2_model'] = 'O2-PWR93'  # 'O2-MPM93'
        # 'N2-SelfContMPM93'
        retrieval_param['n2_model'] = 'N2-SelfContStandardType'

        retrieval_param['selected_species'] = ['O3', retrieval_param['water_vapor_model'],
                                               retrieval_param['o2_model'], retrieval_param['n2_model']]

        #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
        # retrieval_param["azimuth_angle"]=32

        retrieval_param['cira86_path'] = os.path.join(
            ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')
        # Check the structure of the file and maybe use it ?
        # print(netCDF4.Dataset(filename+".nc").groups.keys())

   #     Baseline retrievals
        retrieval_param['sinefit_periods'] = np.array([319e6]) #np.array([319e6])



        assert instrument.instrument_name == instrument_name, 'Wrong instrument definition'

        if recheck_channels:
            integrated_data = instrument.find_bad_channels_stdTb(
                spectrometers=instrument.spectrometers,
                stdTb_threshold=12,
                apply_on='int',
                dimension=['time', 'channel_idx']
            )

        if instrument_name == 'mopi5':
            #instrument = instrument.correction_function_mopi5('U5303', 290)
            integrated_dataset = instrument.correct_troposphere(
                instrument.spectrometers, dim='time', method='Ingold_v1', basis_spectro='AC240')
            instrument.compare_spectra_mopi5(dim='time', idx=[
                                             0, 1, 2, 3], save_plot=False, identifier=[0, 1, 2, 3], with_corr=True)

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
                    retrieval_param['SB_bias'] = 0.05 # in mm
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
                    save_str = 'sensitivity_test_spectroscopy'
                    #line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
                    #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"
                    #line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all_species.par'
                    line_file = '/home/es19m597/Documents/GROMORA/InputsRetrievals/Hitran_f_modified.xml'
                    retrieval_param['line_file'] = line_file
                else: 
                    save_str='other'

                ac, retrieval_param, sensor_out = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, f_bin=None, tb_bin=None, sensor=None)

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
    instrument_name = ['GROMOS','SOMORA']
    date = datetime.date(2018, 6, 9)
    #date = datetime.date(2018, 2, 26)
    cycle = np.arange(7, 8)
    #cycle = np.arange(9, 10)

    tests = ['og', 'noise'] #'continuum','angle','spectroscopy','Tprofile','SB','Tcold','tWindow']
    for radiometer in instrument_name:
        sensitivity_analysis(radiometer, date, param=tests, cycles=cycle)