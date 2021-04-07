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
import time

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from utils_GROSOM import save_single_pdf

from dotenv import load_dotenv

# For ARTS, we need to specify some paths
load_dotenv('/home/es19m597/Documents/ARTS/.env.birg-arts24')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

import mopi5_classes as mc

def retrieve_day(date, instrument_name):
    start = time.time()
   # instrument_name = "mopi5"
    #date = datetime.date(2019,2,22)
    #date = datetime.date(2019,6,15)
   # date = datetime.date(2019,2,16)
    #date = datetime.date(2019,2,1)
    #date = datetime.date(2019,2,6)
    #date = datetime.date(2019,2,13)
    #date = datetime.date(2019,1,4)
    int_time = 6
    integration_strategy = 'classic'
    recheck_channels = True
    do_corr = False
   # cycles = np.arange(0,4)
   # cycles = [8]

    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"

    basename_lvl1 = "/storage/tub/instruments/mopi5/level1/"
    basename_lvl2 = "/home/es19m597/Documents/MOPI5/Data/"

    instrument = mc.MOPI5_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
    
    retrieval_param = instrument.import_standard_retrieval_params_mopi5()

    cycles = np.arange(2,3)

    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param["retrieval_type"] = 2
    retrieval_param['FM_only'] = False
    retrieval_param['show_FM'] = True
    retrieval_param['sensor'] = True
    retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'

    retrieval_param["obs_freq"] = instrument.observation_frequency
    
    retrieval_param["plot_meteo_ds"] = False

    retrieval_param["show_f_grid"] = False
    retrieval_param['plot_opacities'] = False
    retrieval_param["f_shift"] =  0#+1500e3
    retrieval_param["number_of_freq_points"] = 1401
    retrieval_param["irregularity_f_grid"] = 45
    retrieval_param["z_top_sim_grid"] = 112e3
    retrieval_param["z_bottom_sim_grid"] = 600
    retrieval_param["z_resolution_sim_grid"] = 2e3

    retrieval_param["retrieval_grid_type"] = 'altitude'
    retrieval_param["z_top_ret_grid"] = 95e3
    retrieval_param["z_bottom_ret_grid"] = 1e3
    retrieval_param["z_resolution_ret_grid"] = 3e3

    retrieval_param["z_top_ret_grid_h2o"] = 40e3 
    retrieval_param["z_bottom_ret_grid_h2o"] = 600
    retrieval_param["z_resolution_ret_grid_h2o"] = 2e3
    retrieval_param['increased_var_factor'] = 15

    #retrieval_param['unit_var_y']  = 3**2

    retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/es19m597/Documents/GROMORA/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
    retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/es19m597/Documents/GROMORA/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'

    #retrieval_param['obs_freq'] = 1.4217504e11
    retrieval_param['spectroscopy_type'] = 'XML'
    retrieval_param['line_file'] = line_file
    retrieval_param['atm'] ='ecmwf_cira86' # fascod  ecmwf_cira86
    retrieval_param['h2o_apriori']='ecmwf_extended' # 'fascod_extended'
    retrieval_param['ecmwf_store_location'] ='/storage/tub/instruments/gromos/ECMWF_Bern'
    #retrieval_param['ecmwf_store_location'] ='/home/eric/Documents/PhD/ECMWF'
    retrieval_param['extra_time_ecmwf'] = 3.5

    retrieval_param['o3_apriori']='somora'   
    retrieval_param["apriori_O3_cov"] = 1e-6
    retrieval_param["apriori_H2O_stdDev"] = 64e-4 #12e-4

    retrieval_param["apriori_o2_stdDev"]  = 1e-8 #6e-4
    retrieval_param["apriori_n2_stdDev"] = 1e-8

    retrieval_param['water_vapor_model'] = 'H2O-PWR98, H2O' #"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
    retrieval_param['o2_model'] = 'O2-PWR93' #'O2-MPM93'
    retrieval_param['n2_model'] = 'N2-SelfContStandardType' #'N2-SelfContMPM93'
    
    retrieval_param['selected_species']=['O3', retrieval_param['water_vapor_model'],retrieval_param['o2_model'],retrieval_param['n2_model']]
    
    retrieval_param['cira86_path'] = os.path.join(ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')
    retrieval_param['line_file'] = line_file
    # Check the structure of the file and maybe use it ?
    #print(netCDF4.Dataset(filename+".nc").groups.keys())

   # Baseline retrievals 
    retrieval_param['sinefit_periods'] = np.array([319e6])
    retrieval_param['pointing_angle_corr'] = -5
    
    if integration_strategy == 'classic' or integration_strategy == 'TOD_harmo' or integration_strategy == 'TOD':
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(no_flag=False, meta_data=True, extra_base=None)
        dimension=['time', 'channel_idx']
    else:
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(no_flag=True, meta_data=False, extra_base=None)
        dimension=['chunks', 'channel_idx']

    if recheck_channels:
        integrated_data = instrument.find_bad_channels_stdTb(
            spectrometers = instrument.spectrometers, 
            stdTb_threshold = 15,
            apply_on='int',
            dimension=dimension
            )

    if do_corr:
        integrated_dataset = instrument.correct_troposphere(instrument.spectrometers, dim=dimension[0], method='Ingold_v1', basis_spectro='AC240')
    
    # %%
    #
    spectrometers = instrument.spectrometers
  #  spectrometers = ['AC240','U5303']
    spectrometers = ['AC240']
    #spectrometers = ['USRP-A']
   # spectrometers = ['U5303']

    
    #var_fac_U5303_feb = [400, 1500, 1200, 1200, 1500, 1500, 800, 700, 700, 600, 600, 600, 600, 600, 600]
    var_fac_U5303_apr = np.ones((15,1))*100

    #var_factor_USRP_feb = np.ones((15,1))*200
    var_fac_USRP_apr = np.ones((11,1))*20

    var_fac_AC240_apr = np.ones((11,1))*200
    #var_fac_AC240_feb = np.ones((15,1))*600
    #var_fac_AC240[0] = 200

    
    
    #var_factor = {'U5303':var_fac_U5303, 'AC240':var_fac_AC240_feb, 'USRP-A':var_factor_USRP_feb}
    var_factor = {'U5303':var_fac_U5303_apr, 'AC240':var_fac_AC240_apr, 'USRP-A':var_fac_USRP_apr}
    #[450 , 450, 350]

    level2 = dict()
    for i, spectro in enumerate(spectrometers):
        print('Retrieving ', spectro, ' on day: ')
        spectro_dataset = instrument.integrated_data[spectro]

        if spectro == 'USRP-A':
            retrieval_param["number_of_freq_points"] = 2401
            retrieval_param["irregularity_f_grid"] = 90
        figure_list = []
        
        #print(retrieval_param['increased_var_factor'])
        counter = 0
        for c in cycles:
            try:
                retrieval_param['increased_var_factor'] = var_factor[spectro][c]
                retrieval_param["integration_cycle"] = c
           #     instrument.compare_spectra_mopi5(dim=dimension[0], idx=[c], save_plot=False, identifier=np.arange(c+1), with_corr=True)
                #bin_vector = instrument.create_binning(
                #    instrument.frequencies,
                #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
                #    retrieval_param
                #)
                #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)
                if retrieval_param["retrieval_type"] == 1:
                    retrieval_param["surface_altitude"] = 10e3
                    retrieval_param["observation_altitude"] =  15e3
                    ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                    figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(
                        ac, 
                        spectro_dataset, 
                        retrieval_param, 
                        title = 'MOPI5 '+spectro + ' - ' + date.strftime('%Y%m%d')+ ' - ' + str(c), 
                        figure_list = figure_list,
                        fshift=True,
                        bl=True
                    )
                    level2_cycle = ac.get_level2_xarray()
                    #ds1=level2_cycle.expand_dims({'chunks':1})
                    #level2_cycle=level2_cycle.assign_coords({'chunks':c})
                    #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
                    #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
                elif retrieval_param["retrieval_type"] == 2:
                    retrieval_param["surface_altitude"] = 1000
                    retrieval_param["observation_altitude"] =  1000   
                    ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                    if ac.oem_converged:
                        figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='retrieval_o3',figure_list=figure_list)
                        level2_cycle = ac.get_level2_xarray()
                    else:
                        level2_cycle=xr.Dataset()
                    save_str ='_all_retrieved.pdf'
                else:
                    level2_cycle=xr.Dataset()
            except:
                print('problem retrieving cycle : ',c)
            
            if (counter == 0) & (len(level2_cycle) > 0):
                counter = counter + 1
                level2[spectro] = level2_cycle
            elif (counter>0) & (len(level2_cycle) > 0):
                counter = counter + 1
                level2[spectro] = xr.concat([level2[spectro], level2_cycle], dim='time')

        end_time = time.time()
        print("Total time for this day: %s seconds" % (end_time - start))

        if retrieval_param["retrieval_type"] == 1 or retrieval_param["retrieval_type"] == 2:
          #  save_single_pdf(instrument.filename_level2[spectro]+'_'+save_str, figure_list)

            save_single_pdf(instrument.filename_level2[spectro]+'_Perrin_corrected_'+retrieval_param['atm']+'polyfit2.pdf', figure_list)

            #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
            level2[spectro].to_netcdf(path = instrument.filename_level2[spectro]+'.nc')
            # read and plot from the level 2
            #lvl2[spectro] = level2
    
if __name__ == "__main__":
<<<<<<< HEAD
    dates = pd.date_range(start='2019-02-14', end='2019-02-22')
=======
    dates = pd.date_range(start='2019-02-13', end='2019-02-13')
>>>>>>> a23901c5078644d98f5f82093adbf068e47c42c7

    for d in dates:
        try:
            retrieve_day(d, 'mopi5')
        except:
            print('problem retrieving day : ',d)

# %%
