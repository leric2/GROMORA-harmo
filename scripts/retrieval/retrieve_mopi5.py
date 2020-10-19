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

#from dotenv import load_dotenv

import mopi5_classes as mc

# For ARTS, we need to specify some paths
#load_dotenv('/home/eric/Documents/PhD/ARTS/arts-examples/.env.t490-arts2.3')
ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

if __name__ == "__main__":
    instrument_name = "mopi5"
    #date = datetime.date(2019,2,22)
    #date = datetime.date(2019,6,15)
    date = datetime.date(2019,4,27)
    #date = datetime.date(2019,2,1)
    #date = datetime.date(2019,2,6)
    #date = datetime.date(2019,2,13)
    #date = datetime.date(2019,1,4)
    int_time = 1
    integration_strategy = 'meanTb_harmo'
    recheck_channels = True
    do_corr = False
    cycles = np.arange(0,4)

    line_file = ARTS_DATA_PATH+"/spectroscopy/Perrin_newformat_speciessplit/O3-666.xml.gz"

    basename_lvl1 = "/scratch/MOPI5/Level1/"
    basename_lvl2 = "/scratch/MOPI5/Level2/"

    instrument = mc.MOPI5_LvL2(
        date=date,
        basename_lvl1=basename_lvl1,
        basename_lvl2=basename_lvl2,
        integration_strategy=integration_strategy,
        integration_time=int_time
    )
    
    retrieval_param = instrument.import_standard_retrieval_params_mopi5()
    retrieval_param["apriori_O3_cov"] = 1.5e-6
    
    retrieval_param['cira86_path'] = os.path.join(ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly')
    retrieval_param['line_file'] = line_file
    # Check the structure of the file and maybe use it ?
    #print(netCDF4.Dataset(filename+".nc").groups.keys())
    
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
    spectrometers = ['AC240','U5303']
    #spectrometers = ['AC240']
    #spectrometers = ['USRP-A']
    cycles = np.arange(0,11)
    #cycles = [8]
    
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

    retrieval_param['atm'] ='fascod_gromos_o3'

    #lvl2 = dict()
    for i, spectro in enumerate(spectrometers):
        print('Retrieving ', spectro, ' on day: ')
        spectro_dataset = instrument.integrated_data[spectro]

        figure_list = []
        
        #print(retrieval_param['increased_var_factor'])
        counter = 0
        for c in cycles:
            try:
                counter=counter+1
                retrieval_param['increased_var_factor'] = var_factor[spectro][c]
                retrieval_param["integration_cycle"] = c
                instrument.compare_spectra_mopi5(dim=dimension[0], idx=[c], save_plot=False, identifier=np.arange(c+1), with_corr=True)

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
                    retrieval_param["surface_altitude"] = 1200
                    retrieval_param["observation_altitude"] =  1200
                    ac, retrieval_param = instrument.retrieve_cycle(spectro_dataset, retrieval_param, f_bin=None, tb_bin=None)
                    figure_list = instrument.plot_level2(ac, spectro_dataset, retrieval_param, title ='retrieval_o3_h20')
                    save_single_pdf(instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'_Perrin_'+retrieval_param['water_vapor_model']+'.pdf', figure_list)
                else:
                    pass
                
                if counter>1:
                    level2 = xr.concat([level2, level2_cycle], dim='observation')
                else: 
                    level2 = level2_cycle
            except:
                print('could not retrieve this cycle: ',c)
        
        save_single_pdf(instrument.filename_level2[spectro]+'_Perrin_corrected_'+retrieval_param['atm']+'_'+str(cycles)+'.pdf', figure_list)

        #level2.to_netcdf(path = instrument.filename_level2[spectro]+'_'+str(retrieval_param["integration_cycle"])+'.nc')
        level2.to_netcdf(path = instrument.filename_level2[spectro]+'_all_april'+'.nc')
        # read and plot from the level 2
        #lvl2[spectro] = level2

# %%