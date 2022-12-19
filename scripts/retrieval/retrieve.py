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
    instrument_name = "SOMORA"
    date = datetime.date(2021, 1, 20)
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
        import gromos_classes as gc
        basename_lvl1 = os.path.join('/storage/tub/instruments/gromos/level1/GROMORA/v2/', str(date.year))
        instrument = gc.GROMOS_LvL2(
            date,
            basename_lvl1,
            basename_lvl2,
            integration_strategy=integration_strategy,
            integration_time=int_time,
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

    cycles = np.arange(13, 14)

    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'
    retrieval_param['verbose'] = 3
    retrieval_param["retrieval_type"] = 2
    retrieval_param['FM_only'] = False
    retrieval_param['show_FM'] = False
    retrieval_param['date'] = date

    retrieval_param["plot_meteo_ds"] = False

    retrieval_param["show_f_grid"] = False
    retrieval_param['plot_opacities'] = False

    retrieval_param['plot_o3_apriori_covariance'] = True

    retrieval_param = instrument.define_retrieval_param(retrieval_param)

    # Sensor related parameter:
    retrieval_param['sensor'] = 'FFT_SB_Antenna'
    retrieval_param['SB_bias'] = 0
    retrieval_param['FWHM'] = instrument.antenna_fwhm

    # AC240 correction factor
    retrieval_param['AC240_magic_correction'] = True
    retrieval_param["AC240_corr_factor"] = 0.08

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
        elif retrieval_param["retrieval_type"] == 20:
            retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit'
            #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_'+retrieval_param['water_vapor_model']+'.pdf', figure_list)
            retrieval_param["surface_altitude"] = 1000
            retrieval_param["observation_altitude"] = 1000
            retrieval_param["integration_cycle"] = [0, 1, 2]
            ac, retrieval_param = instrument.retrieve_daily(
                spectro_dataset, retrieval_param)
            if ac.oem_converged:
                figure_list = instrument.plot_level2(
                    ac, spectro_dataset, retrieval_param, title='retrieval_o3', figure_list=figure_list)
                level2_cycle = ac.get_level2_xarray()
            else:
                level2_cycle = xr.Dataset()
            save_str = str(
                retrieval_param["integration_cycle"])+'_daily_retrieved.pdf'
        elif retrieval_param["retrieval_type"] == 3:
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] = 1000
            retrieval_param['atm'] = 'ecmwf_cira86'
            retrieval_param['o3_apriori'] = 'gromos'
            retrieval_param['sensor'] = 'FFT_SB_Antenna'
            # retrieval_param['atm']='fascod_gromos_o3'
            retrieval_param['FM_only'] = True
            ac_sim_FM, retrieval_param, sensor_out  = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
            y_noisy = ac_sim_FM.ws.y.value + np.random.normal(0, 0.2, len(ac_sim_FM.ws.f_backend.value))
            retrieval_param['FM_only'] = False
            retrieval_param['sensor'] = 'FFT_SB'
            retrieval_param['o3_apriori'] = 'waccm_monthly'
            ac, retrieval_param, sensor_out  = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, ac_sim_FM=ac_sim_FM, sensor=None, ds_y= y_noisy)
            level2_cycle = ac.get_level2_xarray()
            retrieval_param['sensor'] = 'FFT_SB_Antenna'
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

            ozone_ret, h2o_ret, polyfit_ret, fshift_ret, sinefit_ret  = ac.retrieval_quantities
            ozone_ret2, h2o_ret, polyfit_ret, fshift_ret, sinefit_ret  = ac2.retrieval_quantities
            fig, axs = plt.subplots(nrows=1, ncols=1)
            axs.plot(100*(ozone_ret2.x-ozone_ret.x)/ozone_ret2.x, 1e-3*ozone_ret.z_grid)
            axs.set_xlabel('AT - noAT, [%]')
            axs.set_ylabel('altitude [km]')
            axs.grid()
            fig.savefig('/home/es19m597/Documents/GROMORA/Data/bias_antenna_'+str(date)+'.pdf')

        elif retrieval_param["retrieval_type"] == 4:
            retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] = 800
            # retrieval_param['atm']='fascod_somora_o3'
            # retrieval_param['atm']='fascod_gromos_o3'
            retrieval_param['atm'] = 'ecmwf_cira86'
            retrieval_param['h2o_apriori'] = 'fascod'
            retrieval_param['o3_apriori'] = 'gromos'
            retrieval_param['ref_elevation_angle'] = 90
            retrieval_param['FM_only'] = True
            ac_sim_FM, retrieval_param, sensor_out = instrument.retrieve_cycle(
                spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
            retrieval_param['FM_only'] = False
            # retrieval_param['atm']='fascod_gromos_o3'
            # retrieval_param['atm']='fascod_somora_o3'
            retrieval_param['o3_apriori'] = 'waccm_monthly'
            retrieval_param['h2o_apriori'] = 'ecmwf'
            retrieval_param["number_of_freq_points"] = 5201
            retrieval_param["retrieval_h2o_grid_type"] = 'pressure'
            retrieval_param["h2o_pressure"] = [500e2]
            ac, retrieval_param, sensor_out = instrument.retrieve_cycle(
                spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
            level2_cycle = ac.get_level2_xarray()
            import GROMORA_library
            figure_list1 = GROMORA_library.plot_level2_test_retrieval(
                ac, ac_sim_FM, retrieval_param, title='test_retrieval_o3', z_og=ac_sim_FM.ws.z_field.value[:, 0, 0], og_ozone=ac_sim_FM.ws.vmr_field.value[0, :, 0, 0])
            save_single_pdf(instrument.filename_level2[spectro]+'_'+str(
                c)+'h'+'ecmwf_finer_fgrid'+'.pdf', figure_list1)

            # retrieval_param['h2o_apriori']= 'fascod_extended'
            # ac3, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None, ac=ac_sim_FM)
            # level2_cycle = ac3.get_level2_xarray()
            # figure_list3 = GROSOM_library.plot_level2_test_retrieval(ac3, ac_sim_FM, retrieval_param, title ='test_retrieval_o3', z_og=ac_sim_FM.ws.z_field.value[:,0,0], og_ozone=ac_sim_FM.ws.vmr_field.value[0,:,0,0])
            # save_single_pdf(instrument.filename_level2[spectro]+'_'+str(c)+'h'+'_fascod_extended'+'.pdf', figure_list3)

            # ac2, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None, ac=ac_sim_FM)
            # level2_cycle = ac2.get_level2_xarray()
            # figure_list2 = GROSOM_library.plot_level2_test_retrieval(ac2, retrieval_param, title ='test_retrieval_o3', z_og=ac_sim_FM.ws.z_field.value[:,0,0], og_ozone=ac_sim_FM.ws.vmr_field.value[0,:,0,0])
            # save_single_pdf(instrument.filename_level2[spectro]+'_'+str(c)+'h'+'_H2O-ContMPM93'+'.pdf', figure_list2)
        elif retrieval_param["retrieval_type"] == 18:
            # To compare the FB and the FFT bias
            retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit'
            retrieval_param["surface_altitude"] = 1000
            retrieval_param["observation_altitude"] = 1000
            retrieval_param['o3_apriori'] = 'waccm_monthly'
            retrieval_param['FM_only'] = True
            retrieval_param['sensor']='OFF'
            retrieval_param['sensor'] = 'FFT_SB'
            ac_sim_FM, retrieval_param, sensor_out = instrument.retrieve_cycle(
                spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
            # retrieval_param['atm']='fascod_gromos_o3'
            # retrieval_param['atm']='fascod_somora_o3'
            #retrieval_param['o3_apriori'] = 'waccm_monthly_biased'
            retrieval_param['sensor'] = 'FFT_SB_Antenna'
            ac, retrieval_param, sensor_out = instrument.retrieve_cycle(
                spectro_dataset, retrieval_param, ac_sim_FM=None, sensor=None)
            #level2_cycle = ac.get_level2_xarray()
            fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
            axs[0].plot(1e-9*ac.ws.f_backend.value, ac_sim_FM.y[0], label='No Antenna')
            axs[0].plot(1e-9*ac.ws.f_backend.value, ac.y[0], label='Antenna')
            axs[0].set_ylabel('TB [K]')

            bias= ac.y[0]-  ac_sim_FM.y[0] 
            bias_ds = xr.DataArray(
                data= bias,
                dims='frequency',
                coords={'frequency':ac.ws.f_backend.value}
            )
            bias_ds.rename('bias')
            axin1 = axs[1].inset_axes([0.62, 0.4, 0.35, 0.5])
            inset_ds = bias_ds.where(bias_ds.frequency>142.155*1e9, drop=True).where(bias_ds.frequency<142.195*1e9, drop=True)
            axin1.plot(1e-6*(inset_ds.frequency-instrument.observation_frequency), inset_ds)
            axin1.set_xlabel(r'$\Delta$ f [MHz]')
            axs[1].plot(1e-9*ac.ws.f_backend.value, bias)
            axs[1].set_ylabel(r'$\Delta$ TB [K]')
            axs[1].set_xlabel(r'Frequency [GHz]')
            axs[0].set_title(str(date))
            axs[0].legend()
            #fig.savefig('/home/es19m597/Documents/GROMORA/Data/bias_FB_FFT_'+str(date)+'.pdf')

            plt.show()
            # bias_ds = bias_ds.rename('bias')
            # bias_ds.to_netcdf('/storage/tub/instruments/gromos/spectral_bias_FB-FFT_summer.nc')
            # bias_ds.to_dataframe().to_csv('/storage/tub/instruments/gromos/spectral_bias_FB-FFT_summer.csv')
            # import GROMORA_library

            exit()
            # figure_list1 = GROMORA_library.plot_level2_test_retrieval(
            #     ac, ac_sim_FM, retrieval_param, title='test_retrieval_o3', z_og=ac_sim_FM.ws.z_field.value[:, 0, 0], og_ozone=ac_sim_FM.ws.vmr_field.value[0, :, 0, 0])
            # save_single_pdf(instrument.filename_level2[spectro]+'_'+str(
            #     c)+'h'+'ecmwf_finer_fgrid'+'.pdf', figure_list1)
        elif retrieval_param["retrieval_type"] == 5:
            retrieval_param["surface_altitude"] = 1200
            retrieval_param["observation_altitude"] = 15e3
            ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected_pyarts(
                spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            # figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(
            #     ac,
            #     spectro_dataset,
            #     retrieval_param,
            #     title = 'retrieval_trop_corr',
            #     figure_list = figure_list
            #     )
            #level2_cycle = ac.get_level2_xarray()
            #save_str = str(retrieval_param["integration_cycle"])+'_Perrin_corr.pdf'
        elif retrieval_param["retrieval_type"] == 6:
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] = 800
            ac, retrieval_param = instrument.retrieve_double(
                spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            figure_list = instrument.plot_level2(
                ac, spectro_dataset, retrieval_param, title='retrieval_o3', figure_list=figure_list)
            level2_cycle = ac.get_level2_xarray()
            save_str = str(
                retrieval_param["integration_cycle"])+'_all_retrieved.pdf'
        elif retrieval_param["retrieval_type"] == 7:
            retrieval_param["surface_altitude"] = 900
            retrieval_param["observation_altitude"] = 900
            retrieval_param['FM_only'] = True
            retrieval_param['show_FM'] = False
            retrieval_param['sensor'] = False
            #a_priori = ['mls', 'retrieved_perrin', 'retrieved_hitran']
            a_priori = ['retrieved_gromos', 'retrieved_somora']
            a_priori_legend = ['retrieved_Perrin', 'retrieved_HITRAN']
            o3_apriori_file = [
                '/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/gromos_mean_o3_2017-03-23.nc',
                '/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/somora_mean_o3_2017-03-23.nc',
                # '/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/somora_mean_o3_2017-03-23.nc'
            ]
            ds_freq = spectro_dataset.frequencies[c].values[spectro_dataset.good_channels[c].data == 1]
            y = np.zeros((len(a_priori), 1601))
            fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
            ax1 = axs[0].inset_axes([0.1, 0.4, 0.25, 0.5])

            colors = ['blue', 'red', 'green']
            for i, ap in enumerate(a_priori):
                retrieval_param['o3_apriori'] = ap
                retrieval_param['o3_apriori_file'] = o3_apriori_file[i]
                ac, retrieval_param = instrument.retrieve_cycle(
                    spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                y[i] = ac.y[0]
                ds_freq = ac.f_grid
                axs[0].plot(ds_freq/1e9, y[i],
                            label=a_priori_legend[i], color=colors[i], lw=0.7)
                ax1.plot((ds_freq - instrument.observation_frequency) /
                         1e6, y[i], label=ap, color=colors[i], lw=0.5)
            # for i in[0,1,2]:
            #     ax1.plot(ds_freq, y[i], label = ap, color = colors[i], lw=0.5 )

            axs[1].plot(ds_freq/1e9, y[0]-y[1], label='Per-HIT', color='k')
            #axs[1].plot(ds_freq/1e9, y[0]-y[2], label = 'MLS-HIT', color='green')
            #axs[1].plot(ds_req/1e9, y[0]-y[1], label = 'MLS-Per', color='red')
            axs[0].set_ylabel('Tb [K]')
            #axs[1].set_ylim(-2.5, 2.5)
            ax1.set_xlim(-10, 10)
            ax1.set_ylim(max(y[0])-8, max(y[0])+1)
            axs[1].set_ylabel(r'$\Delta Tb$ [K]')
            axs[1].set_xlabel('frequency [GHz]')
            for ax in axs:
                ax.legend(loc='upper right')
                ax.grid()
            fig.tight_layout()
            plt.show()
            fig.savefig(basename_lvl2+'comparison_FM_sensorFFT_' +
                        instrument.datestr+'.pdf')

            # save_str = str(retrieval_param["integration_cycle"])+'_all_retrieved.pdf'
            level2_cycle = xr.Dataset()
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