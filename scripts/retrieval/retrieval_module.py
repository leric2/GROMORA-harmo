#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:58:02 2020

@author: eric

Module with function for the retrieval of Ozone from GROSOM using pyretrievals

"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

import apriori_data_GROSOM
import GROSOM_library

from retrievals import arts
from retrievals import covmat
from retrievals import utils

from retrievals.arts.atmosphere import p2z_simple, z2p_simple
from retrievals.data import p_interpolate

from typhon.arts.workspace import arts_agenda
#from pyarts.workspace import arts_agenda

def make_f_grid(retrieval_param): #TODO
    '''
    create simulation frequency grid

    '''
    n_f = retrieval_param["number_of_freq_points"]# Number of points
    bw = 1.3*retrieval_param["bandwidth"] # Bandwidth
    x = np.linspace(-1, 1, n_f)
    f_grid = x ** 3 + x / retrieval_param["irregularity_f_grid"]
    f_grid = f_grid * bw / (max(f_grid) - min(f_grid)) + retrieval_param['obs_freq']
    #f_grid = np.linspace(retrieval_param["f_min"]-10, retrieval_param["f_max"]+10, n_f)

    if retrieval_param["show_f_grid"]:
        fig = plt.figure()
        plt.semilogy(f_grid[1:]/1e9,np.diff(f_grid)/1e3,'.')
        plt.xlim((retrieval_param['obs_freq']-200e6)/1e9,(retrieval_param['obs_freq']+200e6)/1e9)
        #plt.ylim(0,300)
        plt.ylabel(r'$\Delta f$ [kHz]')
        plt.suptitle('Frequency grid spacing')
        plt.show()
    return f_grid

def plot_FM_comparison(ds_freq,y_FM,y_obs):
    fig = plt.figure()
    
    fig.suptitle('Comparison between FM and Observation')
    ax = fig.add_subplot(111)
    ax.plot(ds_freq/1e9,y_obs,'r')
    ax.plot(ds_freq/1e9,y_FM,'b-',linewidth=2)
    #ax.set_ylim((-5,25))
    
    ax.set_xlabel('freq [GHz]')
    ax.set_ylabel('Tb [K]')
    
    ax.legend(('Observed','FM'))
    
    plt.show()
    pass

def retrieve_cycle(instrument, spectro_dataset, retrieval_param, ac_FM=None):
    '''
    Retrieval of a single integration cycle defined in retrieval_param

    Parameters
    ----------
    level1b : TYPE
        DESCRIPTION.
    retrieval_param : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    cycle = retrieval_param["integration_cycle"]
    if ac_FM is None:
        print("Retrieval of Ozone and H20")
        good_channels = spectro_dataset.good_channels[cycle].data == 1
        ds_freq = spectro_dataset.frequencies[cycle].values[good_channels]
        ds_y = spectro_dataset.Tb[cycle].values[good_channels]

        ds_num_of_channel = len(ds_freq)
        #ds_Tb = Tb[cycle].values

        ds_bw = max(ds_freq) - min(ds_freq)

        ds_df = ds_bw/(ds_num_of_channel-1)

        retrieval_param["bandwidth"] = instrument.bandwidth[0]
        # defining simulation grids
        if retrieval_param["binned_ch"]:
            f_grid = ds_freq
        else:
            f_grid = make_f_grid(retrieval_param)

        print('Minimum of the frequency grid spacing [kHz]: ', min(np.diff(f_grid))/1e3 )
    else :
        print("Retrieval of Ozone and H20 providing measurement vector")
        ds_freq = ac_FM.ws.f_backend.value
        ds_y = (ac_FM.ws.y.value + 4*np.random.rand(len(ds_freq))) + 5 + 1e-9*(ds_freq-ds_freq[0])*(0) #Gaussian noise + linear baseline possible
        ds_num_of_channel = len(ds_freq)
        #ds_Tb = Tb[cycle].values

        ds_bw = max(ds_freq) - min(ds_freq)

        ds_df = ds_bw/(ds_num_of_channel-1)

        retrieval_param["bandwidth"] = instrument.bandwidth[0]+100e6
        f_grid = make_f_grid(retrieval_param)
        
    # Iniializing ArtsController object
    ac = arts.ArtsController(verbosity=0, agenda_verbosity=0)
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)
    
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
    ac.set_surface(retrieval_param["surface_altitude"])
    
    #spectroscopy
    # abs_species = ["O3","H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252", "O2-PWR93","N2-SelfContStandardType"]
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
        ac.set_atmosphere(atm, vmr_zeropadding=True)
    else:
        ValueError('atmosphere type not recognized')

    # Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
    ac.apply_hse(100e2, 0.5)  # value taken from GROMOS retrieval
    
    # create an observation: 
    # only one for now, but then the whole day ?   
    obs = arts.Observation(
        za = retrieval_param["zenith_angle"], 
        aa = retrieval_param["azimuth_angle"], 
        lat = retrieval_param["lat"],
        lon = retrieval_param["lon"],
        alt = retrieval_param["observation_altitude"],
        time = retrieval_param["time"]
        )

    ac.set_observations([obs])
    
    # Defining our sensors
    sensor = arts.SensorFFT(ds_freq+retrieval_param["f_shift"], ds_df)
    ac.set_sensor(sensor)
    
    # doing the checks
    ac.checked_calc(negative_vmr_ok=False, bad_partition_functions_ok=False)
    
    # FM + noise --> to retrieve as test !
    ac.ws.iy_aux_vars = ["Optical depth"]
    y_FM = ac.y_calc()

    if retrieval_param['FM_only']:
        return ac, retrieval_param

    plot_FM_comparison(ds_freq,y_FM[0],ds_y)
    print('Optical depth values:')
    print(ac.ws.y_aux.value)
    ac.set_y([ds_y])
    
    # Setup the retrieval

    if retrieval_param["retrieval_grid_type"] == 'pressure':
        p_grid_retrieval = np.logspace(5, -1, 61)
    else:
        z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
        z_top_ret = retrieval_param["z_top_ret_grid"]
        z_res_ret = retrieval_param["z_resolution_ret_grid"]
        z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
        p_grid_retrieval = z2p_simple(z_grid_retrieval)
    
    z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
    z_top_ret = retrieval_param["z_top_ret_grid_h2o"]
    z_res_ret = retrieval_param["z_resolution_ret_grid_h2o"]
    z_grid_retrieval_h2o = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
    p_grid_retrieval_h2o = z2p_simple(z_grid_retrieval_h2o)

    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])

    #sx = covmat.covmat_diagonal_sparse(retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval))
    
    sx_water = covmat.covmat_diagonal_sparse(retrieval_param["apriori_H2O_stdDev"] * np.ones_like(p_grid_retrieval_h2o))
    #sx_water = covmat.covmat_diagonal_sparse((z_grid_retrieval_h2o/z_grid_retrieval_h2o[0])*retrieval_param["apriori_H2O_stdDev"])
        
    sx_o2 = covmat.covmat_diagonal_sparse(retrieval_param["apriori_o2_stdDev"] * np.ones_like(p_grid_retrieval))
    sx_n2 = covmat.covmat_diagonal_sparse(retrieval_param["apriori_n2_stdDev"] * np.ones_like(p_grid_retrieval))
    # sx_water = covmat.covmat_1d_sparse(
    #     grid1 = p_grid_retrieval_h2o,
    #     sigma1 = retrieval_param["apriori_H2O_stdDev"] * np.ones_like(p_grid_retrieval_h2o),
    #     cl1 = 0.2 * np.ones_like(p_grid_retrieval),
    #     fname="exp",
    #     cutoff=0.01,
    # )

    sx = covmat.covmat_1d_sparse(
        grid1 = np.log10(p_grid_retrieval),
        sigma1 = retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval),
        cl1 = 0.3 * np.ones_like(p_grid_retrieval),
        fname="lin",
        cutoff=0.001,
    )
    
    ozone_ret = arts.AbsSpecies(
        species = 'O3',
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx,
        unit = 'vmr'
    )
    h2o_ret = arts.AbsSpecies(
        species = water_vapor_model,
        p_grid = p_grid_retrieval_h2o,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx_water,
        unit = 'vmr'
    )
    
    o2_ret = arts.AbsSpecies(
        species = retrieval_param['o2_model'],
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx_o2,
        unit = 'rel'
    )
    n2_ref = arts.AbsSpecies(
        species = retrieval_param['n2_model'] ,
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx_n2,
        unit = 'rel'
    )
    #polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])
    
    # increase variance for spurious spectra by factor
    #factor = retrieval_param['increased_var_factor']
    
    #y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)
    y_var = retrieval_param['increased_var_factor']*np.square(spectro_dataset.noise_level[cycle].data)* np.ones_like(ds_y)
    print('Measurement variance : ', y_var) 
    #y_var[(level1b_ds.good_channels[cycle].values==0)] = factor*retrieval_param['unit_var_y']
    
    #y_var = retrieval_param['increased_var_factor']*np.square(spectro_dataset.stdTb[cycle].data[good_channels])
    polyfit_ret = arts.Polyfit(
       poly_order=1, covmats=[np.array([[5]]), np.array([[1]])]
    )

    fshift_ret = arts.FreqShift(100e3, df=50e3)

    #ac.define_retrieval(retrieval_quantities=[ozone_ret], y_vars=y_var)
    #ac.define_retrieval(retrieval_quantities=[ozone_ret, h2o_ret], y_vars=y_var)
    if retrieval_param['retrieved_h2o']:
        ac.define_retrieval(retrieval_quantities=[ozone_ret, h2o_ret], y_vars=y_var)
    else:
        ac.define_retrieval(retrieval_quantities=[ozone_ret, h2o_ret, o2_ret, n2_ref, fshift_ret],  y_vars=y_var)
    
    # Let a priori be off by 0.5 ppm (testing purpose)
    #vmr_offset = -0.5e-6
    #ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)
    
    # Run retrieval (parameter taken from MOPI)
    # SOMORA is using 'lm': Levenberg-Marquardt (LM) method
    ac.oem(
        method='gn',
        max_iter=10,
        stop_dx=0.1,
        display_progress=1,
        lm_ga_settings = [1000,2,2,1e5,1,99],
        inversion_iterate_agenda=inversion_iterate_agenda,
      )
    
    if not ac.oem_converged:
        print("OEM did not converge.")
        print("OEM diagnostics: " + str(ac.oem_diagnostics))
        for e in ac.oem_errors:
            print("OEM error: " + e)
            continue
    return ac, retrieval_param

@arts_agenda
def inversion_iterate_agenda(ws):
    """Custom inversion iterate agenda to ignore bad partition functions."""
    ws.Ignore(ws.inversion_iteration_counter)

    # Map x to ARTS' variables
    ws.x2artsAtmAndSurf()
    ws.x2artsSensor()

    # To be safe, rerun some checks
    ws.atmfields_checkedCalc(negative_vmr_ok=True)
    ws.atmgeom_checkedCalc()

    # Calculate yf and Jacobian matching x
    ws.yCalc(y = ws.yf)

    # Add baseline term
    ws.VectorAddVector(ws.yf, ws.y, ws.y_baseline)

    # This method takes cares of some "fixes" that are needed to get the Jacobian
    # right for iterative solutions. No need to call this WSM for linear inversions.
    ws.jacobianAdjustAndTransform()

def retrieve_cycle_tropospheric_corrected(instrument, spectro_dataset, retrieval_param, ac_FM=None):
    '''
    Retrieval of a single integration cycle defined in retrieval_param

    Parameters
    ----------
    level1b : TYPE
        DESCRIPTION.
    retrieval_param : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    print("Tropospheric corrected spectra retrieval")

    cycle = retrieval_param["integration_cycle"]
    if ac_FM is None:
        print("Retrieval of Ozone")
        good_channels = spectro_dataset.good_channels[cycle].data == 1
        ds_freq = spectro_dataset.frequencies[cycle].values[good_channels]
        ds_y = spectro_dataset.Tb_corr[cycle].values[good_channels]

        ds_num_of_channel = len(ds_freq)
        #ds_Tb = Tb[cycle].values

        ds_bw = max(ds_freq) - min(ds_freq)

        ds_df = ds_bw/(ds_num_of_channel-1)

        retrieval_param["bandwidth"] = instrument.bandwidth[0]
        # defining simulation grids
        if retrieval_param["binned_ch"]:
            f_grid = ds_freq
        else:
            f_grid = make_f_grid(retrieval_param)
    else :
        print("Retrieval of Ozone providing measurement vector")
        ds_freq = ac_FM.ws.f_backend.value
        ds_y = (ac_FM.ws.y.value + 3*np.random.rand(len(ds_freq))) + 1e-9*(ds_freq-ds_freq[0])*(0) #Gaussian noise + linear baseline possible
        ds_num_of_channel = len(ds_freq)
        #ds_Tb = Tb[cycle].values

        ds_bw = max(ds_freq) - min(ds_freq)

        ds_df = ds_bw/(ds_num_of_channel-1)

        retrieval_param["bandwidth"] = instrument.bandwidth[0]+100e6
        f_grid = make_f_grid(retrieval_param)


    retrieval_param["zenith_angle"] = retrieval_param['ref_elevation_angle'] - spectro_dataset.mean_sky_elevation_angle.values[cycle]
    retrieval_param["azimuth_angle"] = spectro_dataset.azimuth_angle.values[cycle]
    retrieval_param["time"] = spectro_dataset.time[cycle].values
    retrieval_param["lat"] = spectro_dataset.lat[cycle].values
    retrieval_param["lon"] = spectro_dataset.lon[cycle].values
    retrieval_param['time_start'] = spectro_dataset.first_sky_time[cycle].values
    retrieval_param['time_stop'] = spectro_dataset.last_sky_time[cycle].values
    retrieval_param["f_max"] = max(ds_freq)
    retrieval_param["f_min"] = min(ds_freq)
    retrieval_param["bandwidth"] = 1e9

    ds_num_of_channel = len(ds_freq)
    #ds_Tb = Tb[cycle].values
    
    ds_bw = max(ds_freq) - min(ds_freq)
    
    ds_df = ds_bw/(ds_num_of_channel-1)
        
    # Iniializing ArtsController object
    ac = arts.ArtsController(verbosity=0, agenda_verbosity=0)
    ac.setup(atmosphere_dim=1, iy_unit='PlanckBT', ppath_lmax=-1, stokes_dim=1)

    print('Minimum of the frequency grid spacing [kHz]: ', min(np.diff(f_grid))/1e3 )

    z_bottom = retrieval_param["z_bottom_sim_grid"]
    z_top = retrieval_param["z_top_sim_grid"]
    z_res = retrieval_param["z_resolution_sim_grid"]
    z_grid = np.arange(z_bottom, z_top, z_res)
    p_grid = z2p_simple(z_grid)

    ac.set_grids(f_grid, p_grid)
    
    # altitude for the retrieval
    #ac.set_surface(alt.values[cycle])
    ac.set_surface(retrieval_param["surface_altitude"])
    
    # Spectroscopy
    if retrieval_param['spectroscopy_type'] == 'HITRAN': 
        ac.set_spectroscopy_from_file2(
            abs_lines_file=retrieval_param['line_file'],
            abs_species=retrieval_param['selected_species'],
            format='HITRAN',
            line_shape=["VVH", 750e9],
        )
    else:
        ac.set_spectroscopy_from_file(
            abs_lines_file=retrieval_param['line_file'],
            abs_species=retrieval_param['selected_species'],
            format='Arts',
            line_shape=["VVH", 750e9],
        )

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
        ac.set_atmosphere(atm, vmr_zeropadding=True)
    else:
        ValueError('atmosphere type not recognized')
    
    ac.set_atmosphere(atm, vmr_zeropadding=True)
    
    # Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
    ac.apply_hse(100e2, 0.5)  # value taken from GROMOS retrieval

    # create an observation: 
    # only one for now, but then the whole day ?
    obs = arts.Observation(
        za = retrieval_param["zenith_angle"], 
        aa = retrieval_param["azimuth_angle"], 
        lat = retrieval_param["lat"],
        lon = retrieval_param["lon"],
        alt = retrieval_param["observation_altitude"],
        time = retrieval_param["time"]
        )

    ac.set_observations([obs])
    
    # Defining our sensors
    sensor = arts.SensorFFT(ds_freq, ds_df)
    ac.set_sensor(sensor)
    
    # doing the checks
    ac.checked_calc()

    #y_FM = ac.ws.yCalc()
    y_FM = ac.y_calc()

    if retrieval_param['FM_only']:
        return ac, retrieval_param

    plot_FM_comparison(ds_freq,y_FM[0],ds_y)

    ac.set_y([ds_y])
    
    # Setup the retrieval grid
    #VectorNLogSpace( p_ret_grid, np, 500e2, 0.1 )

    if retrieval_param["retrieval_grid_type"] == 'pressure':
        p_grid_retrieval = np.logspace(5, -1, 81)
    else:
        z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
        z_top_ret = retrieval_param["z_top_ret_grid"]
        z_res_ret = retrieval_param["z_resolution_ret_grid"]
        z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
        p_grid_retrieval = z2p_simple(z_grid_retrieval)

    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])
    
    sx = covmat.covmat_diagonal_sparse(retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval))
    # sx = covmat.covmat_1d_sparse(
    #     grid1 = np.log10(p_grid_retrieval),
    #     sigma1 = retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval),
    #     cl1 = 0.1 * np.ones_like(p_grid_retrieval),
    #     fname="lin",
    #     cutoff=0.001,
    # )
    # The different things we want to retrieve
    fshift_ret = arts.FreqShift(100e3, df=50e3)

    ozone_ret = arts.AbsSpecies(
        species = 'O3',
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx,
        unit = 'vmr'
    )
    
    #y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)
    y_var = retrieval_param['increased_var_factor']*np.square(spectro_dataset.noise_level[cycle].data)* np.ones_like(ds_y)
    print('Measurement variance : ', y_var)
    
    #y_var = retrieval_param['increased_var_factor']*np.square(spectro_dataset.stdTb[cycle].data[good_channels])

    #polyfit_ret = arts.Polyfit(
    #    poly_order=1, covmats=[np.array([[5]]), np.array([[1]])]
    #)

    ac.define_retrieval(retrieval_quantities=[ozone_ret], y_vars=y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)

    # # Trying to define_retrieval()
    # ac.ws.retrievalDefInit()

    # ac.ws.covmat_block = []
    # ac.ws.covmat_inv_block = []


    # arts.boilerplate.clear_sparse(ac.ws, ac.ws.covmat_block)
    # arts.boilerplate.clear_sparse(ac.ws, ac.ws.covmat_inv_block)

    # ozone_ret.apply(ac.ws)

    # # ac.ws.retrievalAddAbsSpecies(
    # #     species = 'O3',
    # #     unit='vmr',
    # #     g1 = p_ret_grid,
    # #     g2 = lat_grid,
    # #     g3 = lon_grid
    # # )

    # # arts.boilerplate.set_variable_by_xml(ac.ws,ac.ws.covmat_block, sx_O3)
    # # ac.ws.covmat_sxAddBlock(block=ac.ws.covmat_block)

    # # Se and its inverse
    # from scipy import sparse
    # covmat_block = sparse.diags(y_var, format='csr')
    # arts.boilerplate.set_variable_by_xml(ac.ws, ac.ws.covmat_block, covmat_block)
    # ac.ws.covmat_seAddBlock(block=ac.ws.covmat_block)

    # covmat_block = sparse.diags(1/y_var, format='csr')
    # arts.boilerplate.set_variable_by_xml(ac.ws, ac.ws.covmat_block, covmat_block)
    # ac.ws.covmat_seAddInverseBlock(block=ac.ws.covmat_block)

    # ac.ws.retrievalDefClose()

    # Let a priori be off by 0.5 ppm (testing purpose)
    #vmr_offset = -0.5e-6
    #ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)
    
    # Run retrieval (parameter taken from MOPI)
    # SOMORA is using 'lm': Levenberg-Marquardt (LM) method
    ac.oem(
        method='gn',
        max_iter=10,
        stop_dx=0.01,
        inversion_iterate_agenda=inversion_iterate_agenda,
    )
    
    if not ac.oem_converged:
        print("OEM did not converge.")
        print("OEM diagnostics: " + str(ac.oem_diagnostics))
        for e in ac.oem_errors:
            print("OEM error: " + e)
            continue
    return ac, retrieval_param

def forward_model_tropospheric_corrected(spectro_dataset, retrieval_param):
    '''
    Retrieval of a single integration cycle defined in retrieval_param

    Parameters
    ----------
    level1b : TYPE
        DESCRIPTION.
    retrieval_param : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    print("Tropospheric corrected spectra retrieval")

    cycle = retrieval_param["integration_cycle"]
    good_channels = spectro_dataset.good_channels[cycle].data == 1
    ds_freq = spectro_dataset.frequencies[cycle].values[good_channels]
    ds_y = spectro_dataset.Tb_corr[cycle].values[good_channels]

    retrieval_param["zenith_angle"] = retrieval_param['ref_elevation_angle'] - spectro_dataset.mean_sky_elevation_angle.values[cycle]
    retrieval_param["azimuth_angle"] = spectro_dataset.azimuth_angle.values[cycle]
    retrieval_param["time"] = spectro_dataset.time[cycle].values
    retrieval_param["lat"] = spectro_dataset.lat[cycle].values
    retrieval_param["lon"] = spectro_dataset.lon[cycle].values
    retrieval_param['time_start'] = spectro_dataset.first_sky_time[cycle].values
    retrieval_param['time_stop'] = spectro_dataset.last_sky_time[cycle].values
    retrieval_param["f_max"] = max(ds_freq)
    retrieval_param["f_min"] = min(ds_freq)
    retrieval_param["bandwidth"] = 1e9

    ds_num_of_channel = len(ds_freq)
    #ds_Tb = Tb[cycle].values
    
    ds_bw = max(ds_freq) - min(ds_freq)
    
    ds_df = ds_bw/(ds_num_of_channel-1)
        
    # Iniializing ArtsController object
    ac = arts.ArtsController(verbosity=1, agenda_verbosity=1)
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)

    
    # defining simulation grids
    f_grid = make_f_grid(retrieval_param)

    print('Minimum of the frequency grid spacing [kHz]: ', min(np.diff(f_grid))/1e3 )

    z_bottom = retrieval_param["z_bottom_sim_grid"]
    z_top = retrieval_param["z_top_sim_grid"]
    z_res = retrieval_param["z_resolution_sim_grid"]
    z_grid = np.arange(z_bottom, z_top, z_res)
    p_grid = z2p_simple(z_grid)

    ac.set_grids(f_grid, p_grid)
    
    # altitude for the retrieval
    #ac.set_surface(alt.values[cycle])
    ac.set_surface(retrieval_param["surface_altitude"])
    
    # Spectroscopy
    ac.set_spectroscopy_from_file(
        abs_lines_file=retrieval_param['line_file'],
        abs_species=["O3", "H2O-PWR98", "O2-PWR93", "N2-SelfContStandardType"],
        format='Arts',
        line_shape=["VVH", 750e9],
    )

    #apriori_data_GROSOM.plot_apriori_cira86(retrieval_param)
    if retrieval_param['atm'] == 'fascod':
        atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param)
    elif retrieval_param['atm'] == 'fascod_2':
        ac.ws.AtmRawRead(basename = "planets/Earth/Fascod/{}/{}".format('midlatitude-winter','midlatitude-winter'))
        ac.ws.AtmFieldsCalc()        
    elif retrieval_param['atm'] == 'fascod_gromos_o3':
        atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param) 
    elif retrieval_param['atm'] =='ecmwf_cira86':
        t1 = pd.to_datetime(retrieval_param['time_start'])
        t2 = pd.to_datetime(retrieval_param['time_stop'])
    
        #ecmwf_store = '/home/eric/Documents/PhD/GROSOM/ECMWF'
        ecmwf_store = retrieval_param['ecmwf_store_location'] 
        cira86_path = retrieval_param['cira86_path']

        atm = apriori_data_GROSOM.get_apriori_atmosphere_fascod_ecmwf_cira86(
            retrieval_param,
            ecmwf_store,
            cira86_path,
            t1,
            t2,
            retrieval_param['extra_time_ecmwf']
        )
    else:
        ValueError('atmosphere type not recognized')
    
    ac.set_atmosphere(atm, vmr_zeropadding=True)

    gromos_o3 = apriori_data_GROSOM.read_o3_apriori_ecmwf_mls_gromosOG(retrieval_param['apriori_ozone_climatology_GROMOS'])

    # Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
    ac.apply_hse(100e2, 0.5)  # value taken from GROMOS retrieval

    # create an observation: 
    # only one for now, but then the whole day ?
    obs = arts.Observation(
        za = retrieval_param["zenith_angle"], 
        aa = retrieval_param["azimuth_angle"], 
        lat = retrieval_param["lat"],
        lon = retrieval_param["lon"],
        alt = retrieval_param["observation_altitude"],
        time = retrieval_param["time"]
        )

    ac.set_observations([obs])
    
    # Defining our sensors
    sensor = arts.SensorFFT(ds_freq, ds_df)
    ac.set_sensor(sensor)
    
    # doing the checks
    ac.checked_calc()

    # FM + noise --> to retrieve as test !
    y_FM = ac.y_calc()
    y_FM_noisy = y_FM[0]+5*np.random.rand(len(y_FM[0]))

    plot_FM_comparison(ds_freq,y_FM_noisy,ds_y)

    return y_FM_noisy, gromos_o3

def forward_model(retrieval_param):
    '''
    Function doing a forward model only. This will be the basis for the future 
    retrieve() function

    Parameters
    ----------
    retrieval_param : dict
        Dict containing all types of relevant informations to make the 
        FM. 

    Returns: 
    -------
    None.

    '''
    # Iniializing ArtsController object
    ac = arts.ArtsController()
    
    # atmosphere dim
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)
    
    # defining simulation grids
    f_grid = make_f_grid(retrieval_param) + retrieval_param['obs_freq']
    p_grid = np.logspace(5, -1, 361)
    
    ac.set_grids(f_grid, p_grid)
    
    # altitude for the retrieval
    ac.set_surface(retrieval_param["altitude"])
    
    # spectroscopy
    ac.set_spectroscopy_from_file(
        abs_lines_file = retrieval_param['line_file'],
        abs_species = ["O3", "H2O-PWR98", "O2-PWR98","N2-SelfContStandardType"],
        format = 'Arts',
        line_shape = ("Voigt_Kuntz6", "VVH", 750e9),
        )
    
    ac.set_atmosphere_fascod('midlatitude-winter')
    
    # create an observation:
    fake_obs = arts.Observation(za = retrieval_param["zenith_angle"], aa = retrieval_param["azimuth_angle"], alt = retrieval_param["observation_altitude"])
    
    ac.set_observations([fake_obs])
    
    # Defining our sensors (starting with a fake one)
    #sensor = arts.SensorFFT(ds_f, ds_resolution)
    
    fake_sensor = arts.sensors.SensorOff()
    ac.set_sensor(fake_sensor)
    
    # doing the checks
    ac.checked_calc()
    
    y = ac.y_calc()
    
    return f_grid, y

if __name__=="__main__":
    pass