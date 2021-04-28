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

if __name__ == "__main__":
    start = time.time()
    instrument_name = "mopi5"
    date = datetime.date(2019,2,22)
    #date = datetime.date(2019,6,15)
   # date = datetime.date(2019,2,16)
    #date = datetime.date(2019,2,1)
    #date = datetime.date(2019,2,6)
    #date = datetime.date(2019,2,13)
    #date = datetime.date(2019,1,5)
    int_time = 1
    integration_strategy = 'meanTb_harmo'
    recheck_channels = True
    do_corr = True
   # cycles = np.arange(0,4)
   # cycles = [8]

    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
    line_file = '/home/es19m597/Documents/GROMORA/InputsRetrievals/Hitran_all.xml'

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

    cycles = np.arange(0,15)
    cycles = np.arange(0,9)

    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param["retrieval_type"] = 2
    retrieval_param['FM_only'] = False
    retrieval_param['show_FM'] = False
    retrieval_param['show_profile'] = False
    retrieval_param['sensor'] = True
    retrieval_param['apply_corr_AC240'] = False
    retrieval_param["alpha"] = 0.08
    retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'

    retrieval_param["obs_freq"] = instrument.observation_frequency
    
    retrieval_param["plot_meteo_ds"] = False

    retrieval_param["show_f_grid"] = False
    retrieval_param['plot_opacities'] = False
    retrieval_param["f_shift"] =  0#+1500e3
    retrieval_param["number_of_freq_points"] = 1801
    retrieval_param["irregularity_f_grid"] = 32
    retrieval_param["z_top_sim_grid"] = 100e3
    retrieval_param["z_bottom_sim_grid"] = 600
    retrieval_param["z_resolution_sim_grid"] = 2e3

    retrieval_param["retrieval_grid_type"] = 'altitude'
    retrieval_param["z_top_ret_grid"] = 95e3
    retrieval_param["z_bottom_ret_grid"] = 1e3
    retrieval_param["z_resolution_ret_grid"] = 3e3

    retrieval_param["z_top_ret_grid_h2o"] = 40e3 
    retrieval_param["z_bottom_ret_grid_h2o"] = 600
    retrieval_param["z_resolution_ret_grid_h2o"] = 2e3

    #retrieval_param['unit_var_y']  = 3**2

    retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/es19m597/Documents/GROMORA/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
    retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/es19m597/Documents/GROMORA/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'

    #retrieval_param['obs_freq'] = 1.4217504e11
    retrieval_param['spectroscopy_type'] = 'XML'
    retrieval_param['line_file'] = line_file
    retrieval_param['atm'] ='fascod' # fascod  ecmwf_cira86 fascod_cira86
    retrieval_param['h2o_apriori']='ecmwf_extended' # 'fascod_extended' ecmwf_extended
    retrieval_param['ecmwf_store_location'] ='/storage/tub/instruments/gromos/ECMWF_Bern'
    #retrieval_param['ecmwf_store_location'] ='/home/eric/Documents/PhD/ECMWF'
    retrieval_param['extra_time_ecmwf'] = 3.5

    retrieval_param['add_noise'] = True
    
    
    # For mopi5 chuncks, we define the date where we want ECMWF data
    retrieval_param['time_start'] = datetime.date(2019,1,1)
    retrieval_param['time_stop'] = datetime.date(2019,5,1)

    retrieval_param['o3_apriori']='gromos'   
    retrieval_param["apriori_O3_cov"] = 0.8e-6
    retrieval_param["apriori_H2O_stdDev"] = 12e-4 #12e-4

    retrieval_param["apriori_o2_stdDev"]  = 1e-8 #6e-4
    retrieval_param["apriori_n2_stdDev"] = 1e-8

    retrieval_param['water_vapor_model'] = 'H2O-PWR98, H2O' #"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
    retrieval_param['o2_model'] = 'O2-PWR93' #'O2-MPM93'
    retrieval_param['n2_model'] = 'N2-SelfContStandardType' #'N2-SelfContMPM93'
    
    retrieval_param['selected_species']=['O3', retrieval_param['water_vapor_model'],retrieval_param['o2_model'],retrieval_param['n2_model']]
    
    retrieval_param['cira86_path'] = os.path.join(ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')

    # Check the structure of the file and maybe use it ?
    #print(netCDF4.Dataset(filename+".nc").groups.keys())

   # Baseline retrievals 
    retrieval_param['sinefit_periods'] = np.array([318e6])
    #retrieval_param['sinefit_periods'] = np.array([200e6,318e6])
    retrieval_param['pointing_angle_corr'] = 0
    
    if integration_strategy == 'classic' or integration_strategy == 'TOD_harmo' or integration_strategy == 'TOD':
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(no_flag=False, meta_data=True, extra_base=None)
        dimension=['time', 'channel_idx']
    else:
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(no_flag=True, meta_data=False, extra_base=None)
        dimension=['chunks', 'channel_idx']

    if recheck_channels:
        integrated_data = instrument.find_bad_channels_stdTb(
            spectrometers = instrument.spectrometers, 
            stdTb_threshold = 10,
            apply_on='int',
            dimension=dimension
            )

    if do_corr:
        integrated_dataset = instrument.correct_troposphere(instrument.spectrometers, dim=dimension[0], method='Ingold_v1', basis_spectro='AC240')
    
    # %%
    #
    spectrometers = instrument.spectrometers
    spectrometers = ['AC240_unbiased','AC240']
    
    # spectrometers = ['USRP-A']
    # spectrometers = ['U5303']
    spectrometers = ['U5303','AC240','USRP-A','AC240_unbiased']
    #spectrometers = ['AC240']

    if retrieval_param["retrieval_type"] == 1:
        var_fac_U5303_feb = np.ones((15,1))*18
        var_fac_U5303_feb[2] = 50
        var_factor_USRP_feb = np.ones((15,1))*4
        var_factor_USRP_feb[13] = 8
        var_factor_USRP_feb[14] = 8
        var_fac_AC240_feb = np.ones((15,1))*14
        var_fac_AC240_feb[2] = 38
        var_fac_AC240_feb[13] = 35
        var_fac_AC240_feb[14] = 35
    else:
        if retrieval_param['add_noise']:
            fix_noise_U5303_jan = np.ones((15,1))*0.35
            fix_noise_U5303_jan[4] = 0.5
            fix_noise_AC240_jan = np.ones((15,1))*0.6
            fix_noise_USRP_jan = np.ones((15,1))*0.25

            var_fac_U5303_feb = np.ones((15,1))*1.5
            var_fac_AC240_feb = np.ones((15,1))*1.5
            var_factor_USRP_feb = np.ones((15,1))*1

            fix_noise_U5303_feb = np.ones((15,1))*0.25
            fix_noise_U5303_feb[2] = 0.38
            fix_noise_U5303_feb[[5,6,8]] = 0.3

            fix_noise_AC240_feb = np.ones((15,1))*0.25
            fix_noise_AC240_feb[[0,4,5,6,7,8,10]] = 0.4
            fix_noise_AC240_feb[2] = 0.55
            fix_noise_AC240_feb[1] = 0.15

            fix_noise_USRP_feb = np.ones((15,1))*0.15
            fix_noise_USRP_feb[[11,12,13]] = 0.1
        else:
            var_fac_U5303_feb = np.ones((15,1))*30

            var_fac_U5303_feb[0] = 15
            var_fac_U5303_feb[1] = 10
            var_fac_U5303_feb[2] = 220
            #var_fac_U5303_feb[3] = 10
            var_fac_U5303_feb[4] = 220
            var_fac_U5303_feb[5] = 220
            var_fac_U5303_feb[6] = 60
            var_fac_U5303_feb[7] = 60
            var_fac_U5303_feb[8] = 100
            var_fac_U5303_feb[11] = 20
            var_fac_U5303_feb[12] = 4
            var_fac_U5303_feb[13] = 4
            #var_fac_U5303_feb[14] = 10

            var_factor_USRP_feb = np.ones((15,1))*1
            var_factor_USRP_feb[0] = 3.2
            var_factor_USRP_feb[1] = 3.2
            var_factor_USRP_feb[2] = 3.2
            var_factor_USRP_feb[3] = 2.2
            var_factor_USRP_feb[4] = 3.3
            var_factor_USRP_feb[5] = 1.4
            var_factor_USRP_feb[7] = 1.2

            var_fac_AC240_feb = np.ones((15,1))*25
            var_fac_AC240_feb[0] = 100
            var_fac_AC240_feb[1] = 25
            var_fac_AC240_feb[2] = 220
            # #var_fac_AC240_feb[3] = 1.15
            var_fac_AC240_feb[4] = 220
            var_fac_AC240_feb[5] = 100
            var_fac_AC240_feb[6] = 80
            var_fac_AC240_feb[7] = 80
            var_fac_AC240_feb[8] = 80
            # var_fac_AC240_feb[9] = 8
            # var_fac_AC240_feb[10] = 10
            var_fac_AC240_feb[11] = 16
            var_fac_AC240_feb[12] = 6
            var_fac_AC240_feb[13] = 6
            #var_fac_AC240_feb[14] = 1.15
            #var_fac_AC240[0] = 200
    
    var_factor = {'U5303':var_fac_U5303_feb, 'AC240':var_fac_AC240_feb, 'USRP-A':var_factor_USRP_feb}
    #var_factor = {'U5303':var_fac_U5303_apr, 'AC240':var_fac_AC240_apr, 'USRP-A':var_fac_USRP_apr}
    
    #fix_noise = {'U5303':0.25, 'AC240':0.3, 'USRP-A':0.15}
    fix_noise = {'U5303':fix_noise_U5303_feb, 'AC240':fix_noise_AC240_feb, 'USRP-A':fix_noise_USRP_feb}
    #fix_noise = {'U5303':fix_noise_U5303_jan, 'AC240':fix_noise_AC240_jan, 'USRP-A':fix_noise_USRP_jan}

    #fix_noise = {'U5303':0.25, 'AC240':0.3, 'USRP-A':0.2}

    #var_factor = {'U5303':9, 'AC240':15, 'USRP-A':1.1}
    #var_factor = {'U5303':0.05, 'AC240':0.05, 'USRP-A':0.002}
    #[450 , 450, 350]
    #cycles=[4]
    level2 = dict()
    add_on = '-test'

    THot = instrument.integrated_data['U5303'].THot.data.mean()
    # AC240 bias:
    theoretical_nonlinearities = np.polyfit([80, (THot+80)/2, THot], [0, -0.2, 0], deg=2)
    fitted_poly_theoretical_nonlinearities = np.poly1d(theoretical_nonlinearities)
     
    for ns, spectro in enumerate(spectrometers):
        print('######################################################################################')
        print('Retrieving ', spectro, ' on cycle: ')
        corr = False
        if spectro == 'AC240_unbiased':
            spectro = 'AC240'
            corr = True

        spectro_dataset = instrument.integrated_data[spectro]
        retrieval_param["spectro"] = spectro
        if spectro == 'USRP-A':
            retrieval_param["number_of_freq_points"] = 4601
            retrieval_param["irregularity_f_grid"] = 60
        figure_list = []
        #cycles=np.where(flags[spectro].calibration_flags.data[:,0]==1)[0]
        # if len(cycles) ==0:
        #     return 0
        #print(retrieval_param['increased_var_factor'])
        counter = 0
        for ns, c in enumerate(cycles):
            print(c)
            print('########################################')
            try:
                if corr:
                    good_channels = spectro_dataset.good_channels[c].data == 1
                    ds_freq = spectro_dataset.frequencies[c].values#[good_channels]
                    ds_y = spectro_dataset.Tb[c].values#[good_channels]
                    corr_Tb = (1/(1-retrieval_param["alpha"]))*(ds_y -retrieval_param["alpha"]*np.mean(ds_y))
                    non_lin_delta = fitted_poly_theoretical_nonlinearities(ds_y)
                    corr_Tb_all = (1/(1-retrieval_param["alpha"]))*(ds_y -retrieval_param["alpha"]*np.mean(ds_y) - non_lin_delta)

                    # plt.plot(ds_freq, ds_y, label='non corr')
                    # plt.plot(ds_freq,corr_Tb, label = 'corr')
                    # plt.plot(ds_freq,corr_Tb_all, label = 'corr all')
                    # plt.legend()
                    ds_y = corr_Tb_all
                    spectro_dataset.Tb[c] = ds_y
                    add_on = 'fascod' + 'unbiased_all'
                    

                if retrieval_param['add_noise']:
                    good_channels = spectro_dataset.good_channels[c].data == 1
                    ds_freq = spectro_dataset.frequencies[c].values#[good_channels]
                    ds_y = spectro_dataset.Tb[c].values#[good_channels]
                    noise_y = ds_y + np.random.normal(0, 0, len(ds_freq))    
                    retrieval_param['fix_noise'] = fix_noise[spectro][c]
                    # plt.plot(ds_freq[good_channels],noise_y[good_channels], label = 'corr')
                    # plt.plot(ds_freq[good_channels], ds_y[good_channels], label='non corr')
                    # plt.legend()
                    ds_y = noise_y
                    spectro_dataset.Tb[c] = ds_y

                retrieval_param['increased_var_factor'] = var_factor[spectro][c] #var_factor[spectro][c]

                retrieval_param["integration_cycle"] = c
           #     instrument.compare_spectra_mopi5(dim=dimension[0], idx=[c], save_plot=False, identifier=np.arange(c+1), with_corr=True)
                #bin_vector = instrument.create_binning(
                #    instrument.frequencies,
                #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
                #    retrieval_param
                #)
                #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)
                if retrieval_param["retrieval_type"] == 1:
                    retrieval_param["surface_altitude"] = 1000
                    retrieval_param["observation_altitude"] =  15e3
                    retrieval_param["z_bottom_ret_grid"] = 12e3
                    retrieval_param["z_top_sim_grid"] = 112e3
                    retrieval_param['retrieval_quantities'] = 'o3_fshift_polyfit_sinefit'
                    ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                    # figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(
                    #     ac, 
                    #     spectro_dataset, 
                    #     retrieval_param, 
                    #     title = 'MOPI5 '+spectro + ' - ' + date.strftime('%Y%m%d')+ ' - ' + str(c), 
                    #     figure_list = figure_list,
                    #     fshift=True,
                    #     bl=True
                    # )
                    if ac.oem_converged:
                        figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(
                            ac, 
                            spectro_dataset, 
                            retrieval_param, 
                            title ='retrieval_o3',
                            figure_list=figure_list
                        )
                        level2_cycle = ac.get_level2_xarray()

                    else:
                        level2_cycle=xr.Dataset()
                    save_str ='_all_retrieved.pdf'
                    level2_cycle = ac.get_level2_xarray()
                    #ds1=level2_cycle.expand_dims({'chunks':1})
                    #level2_cycle=level2_cycle.assign_coords({'chunks':c})
                    #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
                    #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_corrected.pdf', figure_list)
                elif retrieval_param["retrieval_type"] == 2:
                    retrieval_param["surface_altitude"] = 900
                    retrieval_param["observation_altitude"] =  900   
                    ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                    if ac.oem_converged:
                        figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='retrieval_o3',figure_list=figure_list)
                        level2_cycle = ac.get_level2_xarray()
                        # yvars = ac.ws.covmat_se.value.to_dense()
                        # var[counter] = yvars[0,0]
                        # level2_cycle['noise_level']= (('time'), )
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

        if retrieval_param["retrieval_type"] == 1:
            save_single_pdf(instrument.filename_level2[spectro]+'_Perrin_corrected_'+retrieval_param['atm']+'polyfit2.pdf', figure_list)
            level2[spectro].to_netcdf(path = instrument.filename_level2[spectro]+'corr.nc')
        else:
            save_single_pdf(instrument.filename_level2[spectro]+'_Perrin_'+retrieval_param['atm']+'polyfit2.pdf', figure_list)
            level2[spectro].to_netcdf(path = instrument.filename_level2[spectro]+add_on+'.nc')
    if retrieval_param['show_profile']:
        plt.plot(ac.ws.t_field.value[:,0,0], ac.ws.z_field.value[:,0,0])
        plt.plot(ac.ws.vmr_field.value[0,:,0,0], ac.ws.z_field.value[:,0,0])
        plt.plot(ac.ws.vmr_field.value[1,:,0,0], ac.ws.z_field.value[:,0,0])
        plt.plot(ac.ws.vmr_field.value[2,:,0,0], ac.ws.z_field.value[:,0,0])
        plt.plot(ac.ws.vmr_field.value[3,:,0,0], ac.ws.z_field.value[:,0,0])

# if __name__ == "__main__":
#     dates = pd.date_range(start='2019-02-22', end='2019-02-22')
#     # date = pd.date_range(start='2019-01-03', end='2019-01-05')
#     # date = pd.date_range(start='2019-01-30', end='2019-02-22')
#     # date = pd.date_range(start='2019-05-01', end='2019-05-04')
#     # No U5303

#     # date = pd.date_range(start='2019-04-25', end='2019-04-27')

#     # date = pd.date_range(start='2019-06-11', end='2019-06-15')

#     for d in dates:
#         try:
#             retrieve_day(d, 'mopi5')
#         except:
#             print('problem retrieving day : ',d)

# %%

