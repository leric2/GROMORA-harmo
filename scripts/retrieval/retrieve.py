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

from dotenv import load_dotenv

# For ARTS, we need to specify some paths
#load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.5')
load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.4')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

if __name__ == "__main__":
    instrument_name = "GROMOS"
    date = datetime.date(2019,1,1)
    int_time = 1
    integration_strategy = 'classic'
    recheck_channels = True
    basename_lvl1 = "/home/eric/Documents/PhD/GROSOM/Data/"
    basename_lvl2 = "/home/eric/Documents/PhD/GROSOM/Data/"
    
    #basename_lvl1 = "/scratch/GROSOM/Level1/"
    #basename_lvl2 = "/scratch/GROSOM/Level2/"
    
    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"
    #line_file = ARTS_DATA_PATH+"/spectroscopy/Hitran/O3-666.xml.gz"
    #line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all_species.par'
    line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all.xml'

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
    
    cycles = np.arange(13,14)
    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param["retrieval_type"] = 2
    retrieval_param['FM_only'] = False
    retrieval_param['retrieved_h2o'] = False
    retrieval_param["obs_freq"] = instrument.observation_frequency
    
    retrieval_param["plot_meteo_ds"] = True
    retrieval_param["number_of_freq_points"] = 1401
    retrieval_param["irregularity_f_grid"] = 45
    retrieval_param["show_f_grid"] = True
    retrieval_param["f_shift"] =  0#+1500e3

    retrieval_param["z_top_sim_grid"] = 97e3
    retrieval_param["z_bottom_sim_grid"] = 600
    retrieval_param["z_resolution_sim_grid"] = 1e3

    retrieval_param["retrieval_grid_type"] = 'altitude'
    retrieval_param["z_top_ret_grid"] = 95e3
    retrieval_param["z_bottom_ret_grid"] = 600
    retrieval_param["z_resolution_ret_grid"] = 2e3

    retrieval_param["z_top_ret_grid_h2o"] = 50e3
    retrieval_param["z_resolution_ret_grid_h2o"] = 1e3

    retrieval_param['increased_var_factor'] = 20

    #retrieval_param['unit_var_y']  = 3**2

    #retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
    #retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'

    #retrieval_param['obs_freq'] = 1.4217504e11
    retrieval_param['spectroscopy_type'] = 'XML'
    retrieval_param['line_file'] = line_file
    retrieval_param['atm'] ='ecmwf_cira86' # fascod  ecmwf_cira86
    retrieval_param['h2o_apriori']='fascod_extended'
    #retrieval_param['ecmwf_store_location'] ='/scratch/ECMWF'
    retrieval_param['ecmwf_store_location'] ='/home/eric/Documents/PhD/ECMWF'
    retrieval_param['extra_time_ecmwf'] = 6

    retrieval_param['o3_apriori']='gromos'   
    retrieval_param["apriori_O3_cov"] = 1e-6
    retrieval_param["apriori_H2O_stdDev"] = 8e-4 #6e-4

    retrieval_param["apriori_o2_stdDev"]  = 1e-4 #6e-4
    retrieval_param["apriori_n2_stdDev"] = 1e-4

    retrieval_param['water_vapor_model'] = 'H2O-MPM93'#"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
    retrieval_param['o2_model'] = 'O2-MPM93' #'O2-MPM93'
    retrieval_param['n2_model'] = 'N2-SelfContMPM93' #'N2-SelfContMPM93'
    
    retrieval_param['selected_species']=['O3',retrieval_param['water_vapor_model'],retrieval_param['o2_model'],retrieval_param['n2_model']]

    #retrieval_param['water_vapor_model'] = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252"
    #retrieval_param["azimuth_angle"]=32

    retrieval_param['cira86_path'] = os.path.join(
        ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')
    # Check the structure of the file and maybe use it ?
    #print(netCDF4.Dataset(filename+".nc").groups.keys())
    
    if integration_strategy == 'classic':
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(no_flag=False, meta_data=True, extra_base=None)
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


    if instrument_name == 'mopi5':
        #instrument = instrument.correction_function_mopi5('U5303', 290)
        integrated_dataset = instrument.correct_troposphere(instrument.spectrometers, dim='time', method='Ingold_v1', basis_spectro='AC240')
        instrument.compare_spectra_mopi5(dim='time', idx=[0,1,2,3], save_plot=False, identifier=[0,1,2,3], with_corr=True)
         
    #for i,s in enumerate(instrument.spectrometer):
    spectro = 'AC240'
    spectro_dataset = instrument.integrated_data[spectro]
    #retrieval_param = {**global_attrs_level1b, **retrieval_param}
    #else :
    #    raise ValueError('incoherent instrument definition')
    
    #level1b_dataset = instrument.smooth_and_apply_correction(level1b_dataset, meteo_ds)
    #instrument.return_bad_channels(date,'U5303')
    #spectro_dataset = instrument.find_bad_channels_stdTb(spectro, 20)

    # instrument.data = instrument.data.find_bad_channels(
    #     spectrometer=spectro,
    #     Tb_min=0,
    #     Tb_max=260,
    #     boxcar_size=128,
    #     boxcar_thresh=7
    # )
    
    # instrument.plot_meteo_ds_level1b_dataset()

    #level1b_dataset = instrument.smooth_corr_spectra(level1b_dataset, retrieval_param)
    #f_sim, y_sim = instrument.forward_model(retrieval_param)
    #plt.plot(f_sim, y_sim[0], level1b_dataset.frequencies.values, level1b_dataset.Tb[1].values)
    

    figure_list = []

    #bin_vector = instrument.create_binning(
    #    instrument.frequencies,
    #    instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data,
    #    retrieval_param
    #)

    #f_bin, tb_bin = instrument.bin_spectrum(instrument.frequencies, instrument.level1b_ds.Tb[retrieval_param["integration_cycle"]].data, bin_vector)
    counter = 0
    for c in cycles:
        counter = counter + 1
        retrieval_param["integration_cycle"] = c
    
        if retrieval_param["retrieval_type"] == 1:
            retrieval_param["surface_altitude"] = 10e3
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
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] =  800
            retrieval_param['plot_opacities'] = True
            ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='retrieval_o3',figure_list=figure_list)
            level2_cycle = ac.get_level2_xarray()
            save_str = str(retrieval_param["integration_cycle"])+'poly_bl.pdf'
            #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_'+retrieval_param['water_vapor_model']+'.pdf', figure_list)
        elif retrieval_param["retrieval_type"] == 3:
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] = 15e3
            retrieval_param['atm']='ecmwf_cira86'
            retrieval_param['o3_apriori']='gromos'
            #retrieval_param['atm']='fascod_gromos_o3'
            retrieval_param['FM_only'] = True
            ac_FM, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            retrieval_param['FM_only'] = False
            retrieval_param['o3_apriori']='somora'
            ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None, ac=ac_FM)
            level2_cycle = ac.get_level2_xarray()
            import GROSOM_library
            figure_list = GROSOM_library.plot_level2_test_retrieval(ac, retrieval_param, title ='test_retrieval_o3', z_og=ac_FM.ws.z_field.value[:,0,0], og_ozone=ac_FM.ws.vmr_field.value[0,:,0,0])
            #save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'Perrin_with_h2o.pdf', figure_list)
        elif retrieval_param["retrieval_type"] == 4:
            retrieval_param["surface_altitude"] = 800
            retrieval_param["observation_altitude"] =  800
            #retrieval_param['atm']='fascod_somora_o3'
            #retrieval_param['atm']='fascod_gromos_o3'
            retrieval_param['atm']='ecmwf_cira86'
            retrieval_param['o3_apriori']='somora'
            retrieval_param['ref_elevation_angle']=90
            retrieval_param['FM_only'] = True
            ac_FM, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            retrieval_param['FM_only'] = False
            #retrieval_param['atm']='fascod_gromos_o3'
            #retrieval_param['atm']='fascod_somora_o3'
            retrieval_param['o3_apriori']='gromos'
            ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None, ac =ac_FM)
            level2_cycle = ac.get_level2_xarray()
            import GROSOM_library
            figure_list = GROSOM_library.plot_level2_test_retrieval(ac, retrieval_param, title ='test_retrieval_o3', z_og=ac_FM.ws.z_field.value[:,0,0], og_ozone=ac_FM.ws.vmr_field.value[0,:,0,0])
            save_single_pdf(instrument.filename_level2[spectro]+'_test'+retrieval_param['water_vapor_model']+'somora_og'+'.pdf', figure_list)
        elif retrieval_param["retrieval_type"] == 5:
            retrieval_param["surface_altitude"] = 1200
            retrieval_param["observation_altitude"] =  15e3
            ac, retrieval_param = instrument.retrieve_cycle_tropospheric_corrected_pyarts(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
            # figure_list = instrument.plot_level2_from_tropospheric_corrected_spectra(
            #     ac, 
            #     spectro_dataset, 
            #     retrieval_param, 
            #     title = 'retrieval_trop_corr',
            #     figure_list = figure_list
            #     )
            #level2_cycle = ac.get_level2_xarray()
            #save_str = str(retrieval_param["integration_cycle"])+'_Perrin_corr.pdf'
        else:
            level2_cycle=[]

        if counter>1:
            level2 = xr.concat([level2, level2_cycle], dim='observation')
        else: 
            level2 = level2_cycle

    #save_single_pdf(instrument.filename_level2[spectro]+'_'+save_str, figure_list)
    #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
