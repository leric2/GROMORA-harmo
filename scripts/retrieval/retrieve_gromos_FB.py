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
from multiprocessing.sharedctypes import Value
import sys, os
from os.path import dirname, abspath, join

sys.path.append(join(dirname(sys.path[0]),'pyretrievals'))
sys.path.append(join(dirname(sys.path[0]),'retrieval'))

import datetime
import time
from abc import ABC

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv

from gromora_utils import save_single_pdf

# For ARTS, you need to specify some paths in your shell. 
# If this is not already done, you can import them now with dotenv:
#load_dotenv('/opt/anaconda/.env.birg-arts24')
#load_dotenv('/opt/arts/.env.stockhorn-arts24')

ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

if __name__ == "__main__":
    start = time.time()
    instrument_name = "GROMOS"
    date = datetime.date(2006, 1, 30)
    int_time = 1
    integration_strategy = 'classic'
    recheck_channels = False

    basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_pyarts/"
    basename_lvl2 = "/home/es19m597/Documents/GROMORA/Data/"

    # Dictionnary containing all EXTERNAL retrieval parameters
    retrieval_param = dict()
    retrieval_param["GROMORA_FOLDER"] = dirname(dirname(dirname(abspath(__file__))))
    retrieval_param["ARTS_DATA_PATH"] = ARTS_DATA_PATH
    retrieval_param["ARTS_BUILD_PATH"] = ARTS_BUILD_PATH
    retrieval_param["ARTS_INCLUDE_PATH"] = ARTS_INCLUDE_PATH

    if instrument_name == "GROMOS":
        import gromos_FB_classes as fb
        basename_lvl1 = os.path.join(
            '/storage/atmosphere/instruments/gromos/level1/GROMORA/v2/', str(date.year))
        instrument = fb.GROMOS_FB_LvL2(
            date,
            basename_lvl1,
            basename_lvl2,
            integration_strategy,
            int_time,
            extra_base=''
            )
    elif instrument_name == "SOMORA":
        raise ValueError('Only GROMOS had FB !')

    if integration_strategy == 'classic':
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
            no_flag=False, meta_data=True, extra_base='')
    else:
        raise NotImplementedError(
            'TODO, implement reading level1b in non classical cases !')

    cycles = np.arange(13,14)

    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit'
    retrieval_param['verbose'] = 3
    retrieval_param["retrieval_type"] = 2
    retrieval_param['FM_only'] = False
    retrieval_param['show_FM'] = True
    retrieval_param['date'] = date

    retrieval_param["plot_meteo_ds"] = False

    retrieval_param["show_f_grid"] = True
    retrieval_param['plot_opacities'] = False

    retrieval_param['plot_o3_apriori_covariance'] = True

    retrieval_param = instrument.define_retrieval_param(retrieval_param)

    # Sensor related parameter:
    retrieval_param['sensor'] = 'FB_SB_Antenna'
    retrieval_param['SB_bias'] = 0
    retrieval_param['FWHM'] = instrument.antenna_fwhm

    # AC240 correction factor
    retrieval_param['AC240_magic_correction'] = False
    retrieval_param["AC240_corr_factor"] = 0

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
    spectro = 'FB'
    spectro_dataset = instrument.integrated_data[spectro]
    if date <  datetime.date(1994, 10 , 1):
        spectro_dataset['mean_sky_elevation_angle'] = ('time',55*np.ones_like(spectro_dataset.mean_sky_elevation_angle.data))
    
    figure_list = []

    # bin_vector = instrument.create_binning(
    #    instrument.frequencies,
    #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
    #    retrieval_param
    # )

    #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)
    counter = 0
    for c in cycles:
        counter = counter + 1
        retrieval_param["integration_cycle"] = c
        print('retrieving cycle : ', c)

        if ~np.isnan(integrated_meteo[spectro].air_pressure[c].data):
            retrieval_param['p_surface'] = integrated_meteo[spectro].air_pressure[c].data
        else:
            retrieval_param['p_surface'] = instrument.standard_air_pressure
        
        if ~np.isnan(integrated_meteo[spectro].air_temperature[c].data):
            retrieval_param['T_surface'] = integrated_meteo[spectro].air_temperature[c].data
        else:
            retrieval_param['T_surface'] = instrument.standard_air_temperature

        if retrieval_param["retrieval_type"] == 1:
            retrieval_param["surface_altitude"] = 1e3
            retrieval_param["observation_altitude"] = 15e3
            ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(
                spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(
                ac,
                spectro_dataset,
                retrieval_param,
                title='retrieval_trop_corr',
                figure_list=figure_list
            )
            level2_cycle = ac.get_level2_xarray()
            save_str = str(
                retrieval_param["integration_cycle"])+'_Perrin_corr.pdf'
            #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
            #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
        elif retrieval_param["retrieval_type"] == 2:
            retrieval_param["surface_altitude"] = 1000
            retrieval_param["observation_altitude"] = 1000
            if c > 1:
                ac, retrieval_param, sensor_out = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
            else:
                ac, retrieval_param, sensor_out = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, ac_sim_FM=None)
            if ac.oem_converged:
                figure_list = instrument.plot_level2(
                    ac, spectro_dataset, retrieval_param, title='retrieval_o3', figure_list=figure_list)
               # level2_cycle = xr.Dataset()
                #ac.level2_diagnostics()
                level2_cycle = ac.get_level2_xarray()
            else:
                level2_cycle = xr.Dataset()
            save_str = str(
                retrieval_param["integration_cycle"])+'_all_retrieved.pdf'
        elif retrieval_param["retrieval_type"] == 3:
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] = 15e3
            retrieval_param['atm'] = 'ecmwf_cira86'
            retrieval_param['o3_apriori'] = 'gromos'
            # retrieval_param['atm']='fascod_gromos_o3'
            retrieval_param['FM_only'] = True
            ac_sim_FM, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(
                spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            retrieval_param['FM_only'] = False
            retrieval_param['o3_apriori'] = 'somora'
            ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(
                spectro_dataset, retrieval_param, f_bin=None, tb_bin=None, ac=ac_sim_FM)
            level2_cycle = ac.get_level2_xarray()
            import retrieval.GROMORA_library as GROMORA_library
            figure_list = GROMORA_library.plot_level2_test_retrieval(
                ac, retrieval_param, title='test_retrieval_o3', z_og=ac_sim_FM.ws.z_field.value[:, 0, 0], og_ozone=ac_sim_FM.ws.vmr_field.value[0, :, 0, 0])
            #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'Perrin_with_h2o.pdf', figure_list)
        elif retrieval_param["retrieval_type"] == 4:
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] = 1000
            retrieval_param['atm'] = 'era5_cira86'
            retrieval_param['o3_apriori'] = 'gromos'
            retrieval_param['sensor'] = 'FB_SB_Antenna'
            # Frequency grid for the simulation: (only for simulated spectra)
            retrieval_param["number_of_freq_points"] = 1201
            retrieval_param["irregularity_f_grid"] = 45
            # retrieval_param['atm']='fascod_gromos_o3'
            retrieval_param['FM_only'] = True
            ac_sim_FM, retrieval_param, sensor_out  = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
            y_noisy = ac_sim_FM.ws.y.value + np.random.normal(0, 0.2, len(ac_sim_FM.ws.f_backend.value))
            retrieval_param['FM_only'] = False
            retrieval_param['sensor'] = 'FB_SB'
            retrieval_param['o3_apriori'] = 'waccm_monthly'
            ac, retrieval_param, sensor_out  = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, ac_sim_FM=ac_sim_FM, sensor=None, ds_y= y_noisy)
            level2_cycle = ac.get_level2_xarray()
            retrieval_param['sensor'] = 'FB_SB_Antenna'
            ac2, retrieval_param, sensor_out  = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, ac_sim_FM=ac_sim_FM, sensor=None, ds_y= y_noisy)
            import GROMORA_library
            figure_list = GROMORA_library.plot_level2_test_retrieval(
                ac, ac_sim_FM, retrieval_param, title='test_retrieval_o3', z_og=ac_sim_FM.ws.z_field.value[:, 0, 0], og_ozone=ac_sim_FM.ws.vmr_field.value[0, :, 0, 0])
            #plt.show()
            save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'generic_no_antenna.pdf', figure_list)
            figure_list2 = GROMORA_library.plot_level2_test_retrieval(
                ac2, ac_sim_FM, retrieval_param, title='test_retrieval_o3', z_og=ac_sim_FM.ws.z_field.value[:, 0, 0], og_ozone=ac_sim_FM.ws.vmr_field.value[0, :, 0, 0])
            save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'generic_antenna.pdf', figure_list2)

            ozone_ret, h2o_ret, polyfit_ret, fshift_ret  = ac.retrieval_quantities
            ozone_ret2, h2o_ret, polyfit_ret, fshift_ret  = ac2.retrieval_quantities
            fig, axs = plt.subplots(nrows=1, ncols=1)
            axs.plot(100*(ozone_ret2.x-ozone_ret.x)/ozone_ret2.x, 1e-3*ozone_ret.z_grid)
            axs.set_xlabel('AT - noAT, [%]')
            axs.set_ylabel('altitude [km]')
            axs.grid()
            fig.savefig('/home/es19m597/Documents/GROMORA/Data/bias_antenna_'+str(date)+'.pdf')
        else:
            level2_cycle = xr.Dataset()

        if counter > 1:
            level2 = xr.concat([level2, level2_cycle], dim='time')
        else:
            level2 = level2_cycle
        # except:
        #    print('problem retrieving cycle : ',c)
        end_time = time.time()
        print("Total time for this cycle: %s seconds" % (end_time - start))
        start = time.time()

    if retrieval_param["retrieval_type"] == 1 or retrieval_param["retrieval_type"] == 2:
        save_single_pdf(
            instrument.filename_level2[spectro]+'_'+save_str, figure_list)

        level2 = instrument.write_level2_gromora(level2, retrieval_param, full_name = instrument.filename_level2[instrument.spectrometers[0]]+'_'+str(retrieval_param["integration_cycle"])+'.nc')