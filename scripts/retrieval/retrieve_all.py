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

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
import gc

from utils_GROSOM import save_single_pdf

# For ARTS, we need to specify some paths
load_dotenv('/opt/anaconda/.env.birg-arts24_pyarts')


ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

def retrieve_day(date, instrument_name):
  #  instrument_name = "SOMORA"
   # date = datetime.date(2019,4,d)
    int_time = 1
    integration_strategy = 'classic'
    recheck_channels = False

    # Dictionnary containing all EXTERNAL retrieval parameters 
    retrieval_param = dict()
    retrieval_param["ARTS_DATA_PATH"] = ARTS_DATA_PATH
    retrieval_param["ARTS_BUILD_PATH"] = ARTS_BUILD_PATH
    retrieval_param["ARTS_INCLUDE_PATH"] = ARTS_INCLUDE_PATH


    if instrument_name=="GROMOS":
        import gromos_classes as gc
        basename_lvl1 = os.path.join('/storage/tub/instruments/gromos/level1/GROMORA/v2/',str(date.year))
        basename_lvl2 = os.path.join('/storage/tub/instruments/gromos/level2/GROMORA/v2/',str(date.year))
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
        basename_lvl1 = os.path.join('/storage/tub/instruments/somora/level1/v2/',str(date.year))
        basename_lvl2 = os.path.join('/storage/tub/instruments/somora/level2/v2/',str(date.year))
        import somora_classes as sm
        instrument = sm.SOMORA_LvL2(
            date=date,
            basename_lvl1=basename_lvl1,
            basename_lvl2=basename_lvl2,
            integration_strategy=integration_strategy,
            integration_time=int_time,
            extra_base=''
        )
        retrieval_param['increased_var_factor'] = 1 # 1.1 for constant o3 cov 
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
    
    if integration_strategy == 'classic':
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
            no_flag=False, meta_data=True, extra_base='')
    else:
        raise NotImplementedError(
            'TODO, implement reading level1b in non classical cases !')

   # cycles = np.arange(1,24)
    
    #cycles = [17]
    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit_sinefit'
    retrieval_param['verbose'] = 1
    retrieval_param["retrieval_type"] = 2
    retrieval_param['FM_only'] = False
    retrieval_param['show_FM'] = False
    retrieval_param['date'] = date

    retrieval_param["plot_meteo_ds"] = False

    retrieval_param["show_f_grid"] = False
    retrieval_param['plot_opacities'] = False

    retrieval_param['plot_o3_apriori_covariance'] = False

    retrieval_param = instrument.define_retrieval_param(retrieval_param)

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

    cycles=np.where(flags[spectro].calibration_flags.data[:,0]==1)[0]
    #cycles= cycles[np.where(np.mod(cycles,8)==1)[0]]
   # cycles = [1,4,10, 15]
    #cycles = [2, 14]
    #cycles= cycles[np.where(np.mod(cycles,8)==1)[0]]
    if len(cycles) ==0:
        return 0
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
               # save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_'+'.pdf', figure_list)
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
                ac_FM, retrieval_param, sensor = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                retrieval_param['FM_only'] = False
                #retrieval_param['atm']='fascod_gromos_o3'
                #retrieval_param['atm']='fascod_somora_o3'
                retrieval_param['o3_apriori']='gromos'
                ac, retrieval_param, sensor = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None, ac =ac_FM)
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
        gc.collect()
        #return None
    # else:
    #     return 0
if __name__ == "__main__":
    void_date_problem = []
        # datetime.date(2009,11,3), 
        # datetime.date(2009,11,27),
        # datetime.date(2009,12,25),
        # datetime.date(2009,12,26), 
        # datetime.date(2010,1,6), 
        # datetime.date(2010,1,7), 
        # datetime.date(2010,1,13),
        # datetime.date(2010,1,31),
        # datetime.date(2010,2,2),
        # datetime.date(2010,3,3),
        # datetime.date(2010,5,31),
    #     datetime.date(2010,7,5),
    #     datetime.date(2010,7,12),
    #     datetime.date(2010,7,29),
    #     datetime.date(2010,8,12),
    #     datetime.date(2010,9,7),
    #     datetime.date(2010,9,16),
    #     datetime.date(2010,9,18),
    #     datetime.date(2010,9,19)]
        #datetime.date(2017,5,26), 
        # datetime.date(2018,5,5),
        # datetime.date(2018,12,24),
        # datetime.date(2018,12,25), 
        # datetime.date(2018,12,26), 
        #datetime.date(2019,1,3)]
    # dates = [
    #     datetime.date(2014,8,12),
    #     datetime.date(2014,9,2),
    #     datetime.date(2010,12,6),
    #     datetime.date(2010,12,8),
    #     datetime.date(2010,12,9),
    #     datetime.date(2011,6,9),
    #     datetime.date(2011,6,30),
    #     datetime.date(2012,4,2),
    #     datetime.date(2012,5,15),
    #     datetime.date(2012,7,3),
    #     datetime.date(2013,7,30),
    #     datetime.date(2016,10,10),
    #     datetime.date(2017,4,24),
    #     datetime.date(2018,8,23),
    #     datetime.date(2018,10,16),
    #     ]
    dates = [datetime.date(2015,7,6),datetime.date(2015,8,13)]
    #dates = pd.date_range(start=sys.argv[1], end=sys.argv[2])#.append(pd.date_range(start='2014-01-01', end='2014-12-31'))
    #dates = pd.date_range(start='2010-01-06', end='2010-01-07').append(pd.date_range(start='2010-01-13', end='2010-01-13')).append(pd.date_range(start='2010-05-31', end='2010-05-31'))#.append(pd.date_range(start='2012-11-26', end='2012-12-31'))#).append(pd.date_range(start='2016-12-31', end='2017-01-01'))
    print('######################################################################################')
    print('######################################################################################')
    print('######################################################################################')
    for d in dates:
        if d in void_date_problem :
            print('abort core problem with this day : ',d ,' --> skipping')
        else:
            try:
                retrieve_day(d, 'SOMORA')
            except:
                #print('problem retrieving day : ',d)
                pass
            print('######################################################################################')
            # print('######################################################################################')
            # try:
            #     retrieve_day(d, 'GROMOS')
            # except:
            #     pass
            #     print('problem retrieving day : ',d)
            

