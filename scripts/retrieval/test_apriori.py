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

import apriori_data_GROSOM
import GROSOM_library

from retrievals import arts
from retrievals import covmat
from retrievals import utils

from retrievals.arts.atmosphere import p2z_simple, z2p_simple
from retrievals.data import p_interpolate

from typhon.arts.workspace import arts_agenda
#from pyarts.workspace import arts_agenda



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
    #line_file = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/Hitran_all.xml'

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

    # Dictionnary containing all EXTERNAL retrieval parameters 
    retrieval_param = dict()


    # type of retrieval to do:
    # 1. tropospheric corrected
    # 2. with h20
    # 3. test retrieving the FM
    retrieval_param["retrieval_type"] = 2
    retrieval_param['FM_only'] = False
    retrieval_param['retrieval_quantities'] = 'o3_h2o_fshift_polyfit'

    retrieval_param["obs_freq"] = instrument.observation_frequency
    
    retrieval_param["plot_meteo_ds"] = False

    retrieval_param["show_f_grid"] = False
    retrieval_param['plot_opacities'] = False
    retrieval_param["f_shift"] =  0#+1500e3
    retrieval_param["number_of_freq_points"] = 1201
    retrieval_param["irregularity_f_grid"] = 45
    retrieval_param["z_top_sim_grid"] = 112e3
    retrieval_param["z_bottom_sim_grid"] = 600
    retrieval_param["z_resolution_sim_grid"] = 500

    retrieval_param["retrieval_grid_type"] = 'altitude'
    retrieval_param["z_top_ret_grid"] = 95e3
    retrieval_param["z_bottom_ret_grid"] = 1e3
    retrieval_param["z_resolution_ret_grid"] = 2e3

    retrieval_param["z_top_ret_grid_h2o"] = 40e3 
    retrieval_param["z_bottom_ret_grid_h2o"] = 600
    retrieval_param["z_resolution_ret_grid_h2o"] = 1e3
    retrieval_param['increased_var_factor'] = 15

    #retrieval_param['unit_var_y']  = 3**2

    #retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'
    #retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/esauvageat/Documents/GROSOM/Analysis/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    retrieval_param['apriori_ozone_climatology_SOMORA'] = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    retrieval_param['apriori_ozone_climatology_GROMOS'] = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/apriori_ECMWF_MLS.O3.aa'

    #retrieval_param['obs_freq'] = 1.4217504e11
    retrieval_param['spectroscopy_type'] = 'XML'
    retrieval_param['line_file'] = line_file
    retrieval_param['atm'] ='ecmwf_cira86' # fascod  ecmwf_cira86
    retrieval_param['h2o_apriori']='ecmwf_extended' # 'fascod_extended'
    #retrieval_param['ecmwf_store_location'] ='/scratch/ECMWF'
    retrieval_param['ecmwf_store_location'] ='/home/eric/Documents/PhD/ECMWF'
    retrieval_param['extra_time_ecmwf'] = 6

    retrieval_param['o3_apriori']='somora'   
    retrieval_param["apriori_O3_cov"] = 1e-6
    retrieval_param["apriori_H2O_stdDev"] = 12e-4 #6e-4

    retrieval_param["apriori_o2_stdDev"]  = 1e-8 #6e-4
    retrieval_param["apriori_n2_stdDev"] = 1e-8

    retrieval_param['water_vapor_model'] = 'H2O-PWR98, H2O' #"H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252" #'H2O-MPM93'
    retrieval_param['o2_model'] = 'O2-PWR93' #'O2-MPM93'
    retrieval_param['n2_model'] = 'N2-SelfContStandardType' #'N2-SelfContMPM93'
    
    retrieval_param['selected_species']=['O3', retrieval_param['water_vapor_model'],retrieval_param['o2_model'],retrieval_param['n2_model']]

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

    spectro = 'AC240'
    spectro_dataset = instrument.integrated_data[spectro]

    cycle = 7
    retrieval_param['ref_elevation_angle'] = 90
    
    retrieval_param["observation_altitude"]=800

    ac = arts.ArtsController(verbosity=0, agenda_verbosity=0)
    ac.setup(atmosphere_dim=1, iy_unit='PlanckBT', ppath_lmax=-1, stokes_dim=1)


    good_channels = spectro_dataset.good_channels[cycle].data == 1
    ds_freq = spectro_dataset.frequencies[cycle].values[good_channels]
    ds_y = spectro_dataset.Tb[cycle].values[good_channels]
    ds_num_of_channel = len(ds_freq)
    #ds_Tb = Tb[cycle].values
    ds_bw = max(ds_freq) - min(ds_freq)
    ds_df = ds_bw/(ds_num_of_channel-1)
    retrieval_param["bandwidth"] = instrument.bandwidth[0]
    # defining simulation grids
    f_grid = ds_freq
    
    retrieval_param["zenith_angle"] = retrieval_param['ref_elevation_angle'] - spectro_dataset.mean_sky_elevation_angle.values[cycle]
    retrieval_param["azimuth_angle"] = spectro_dataset.azimuth_angle.values[cycle]
    retrieval_param["time"] = spectro_dataset.time[cycle].values
    retrieval_param["lat"] = spectro_dataset.lat[cycle].values
    retrieval_param["lon"] = spectro_dataset.lon[cycle].values
    retrieval_param['time_start'] = spectro_dataset.first_sky_time[cycle].values
    retrieval_param['time_stop'] = spectro_dataset.last_sky_time[cycle].values
    retrieval_param["f_max"] = max(ds_freq)
    retrieval_param["f_min"] = min(ds_freq)
    
    z_bottom = retrieval_param["z_bottom_sim_grid"]
    z_top = retrieval_param["z_top_sim_grid"]
    z_res = retrieval_param["z_resolution_sim_grid"]
    z_grid = np.arange(z_bottom, z_top, z_res)
    p_grid = z2p_simple(z_grid)

    ac.set_grids(f_grid, p_grid)
    
    # altitude for the retrieval
    #ac.set_surface(level1b_ds.alt.values[cycle])

    
    #spectroscopy
    # abs_species = ["O3","H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252", "O2-PWR93","N2-SelfContStandardType"]
    water_vapor_model = retrieval_param['water_vapor_model'] 

    #for i,s in enumerate(instrument.spectrometer):
    spectro = 'AC240'
    spectro_dataset = instrument.integrated_data[spectro]

    water_vapor_model = retrieval_param['water_vapor_model'] 

    ac.set_spectroscopy_from_file(
        abs_lines_file = retrieval_param['line_file'],
        abs_species = retrieval_param['selected_species'],
        format = 'Arts',
        line_shape = ["VVH", 750e9],
        )
    
    #apriori_data_GROSOM.plot_apriori_cira86(retrieval_param)
    if retrieval_param['atm'][0:6] == 'fascod':
        atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param)
        ac.set_atmosphere(atm, vmr_zeropadding=True)
    elif retrieval_param['atm'] == 'fascod_2':
        ac.ws.AtmRawRead(basename = "planets/Earth/Fascod/{}/{}".format('midlatitude-winter','midlatitude-winter'))
        ac.ws.AtmFieldsCalc()
    elif retrieval_param['atm'] =='ecmwf_cira86':
        t1 = pd.to_datetime(retrieval_param['time_start'])
        t2 = pd.to_datetime(retrieval_param['time_stop'])
    
        #ecmwf_store = '/home/eric/Documents/PhD/GROSOM/ECMWF'
        ecmwf_store = retrieval_param['ecmwf_store_location'] 
        cira86_path = retrieval_param['cira86_path']

        ecmwf_prefix = f'ecmwf_oper_v{2}_{instrument.location}_%Y%m%d.nc'
        retrieval_param['ecmwf_prefix']=ecmwf_prefix

        atm = apriori_data_GROSOM.get_apriori_atmosphere_fascod_ecmwf_cira86(
            retrieval_param,
            ecmwf_store,
            cira86_path,
            t1,
            t2,
            retrieval_param['extra_time_ecmwf']
        )

        print('Min altitude before setting atm: ', min(atm.z_field))
        print('Max altitude before setting atm: ', max(atm.z_field))

        try:
            ac.set_atmosphere(atm, vmr_zeropadding=True)
        except:
            if max(atm.z_field) > max(p_grid):
                p_grid = np.delete(p_grid, 0)
                ac.ws.p_grid = p_grid
                ac.set_atmosphere(atm, vmr_zeropadding=True)
    else:
        ValueError('atmosphere type not recognized')

    print('Min altitude : ', min(ac.ws.z_field.value[:,0,0]))
    print('Max altitude : ', max(ac.ws.z_field.value[:,0,0]))


    retrieval_param["surface_altitude"] = 600

    if retrieval_param["surface_altitude"] < min(ac.ws.z_field.value[:,0,0]):
        retrieval_param["surface_altitude"] = min(ac.ws.z_field.value[:,0,0]) 
        print('Surface altitude has been changed to min of z_field : ', retrieval_param["surface_altitude"])
   
    ac.set_surface(retrieval_param["surface_altitude"])

    # Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
    ac.apply_hse(ac.ws.p_grid.value[0], 0.5)# value taken from GROMOS retrieval

    print('Min altitude (HSE): ', min(ac.ws.z_field.value[:,0,0]))
    print('Max altitude (HSE): ', max(ac.ws.z_field.value[:,0,0]))


    obs = arts.Observation(
    za = retrieval_param["zenith_angle"], 
    aa = retrieval_param["azimuth_angle"], 
    lat = retrieval_param["lat"],
    lon = retrieval_param["lon"],
    alt = retrieval_param["observation_altitude"],
    time = retrieval_param["time"]
    )

    ac.set_observations([obs])
    
    # doing the checks
    negative_vmr_ok = 1

    ac.ws.abs_xsec_agenda_checkedCalc()
    ac.ws.propmat_clearsky_agenda_checkedCalc()
    ac.ws.atmfields_checkedCalc(
        negative_vmr_ok=1 if negative_vmr_ok else 0
        )
    ac.ws.atmgeom_checkedCalc()
    ac.ws.cloudbox_checkedCalc()