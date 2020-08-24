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

import retrieval_module
import apriori_data_GROSOM
import data_GROSOM

from retrievals import arts
from retrievals import covmat
from retrievals import utils

from retrievals.arts.atmosphere import p2z_simple, z2p_simple
from retrievals.data import p_interpolate

from typhon.arts.workspace import arts_agenda

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
        plt.suptitle('Simulation f_grid spacing')
        plt.show()
    return f_grid

def plot_FM_comparison(ds_freq,y_FM,y_obs):
    fig = plt.figure()
    
    fig.suptitle('Comparison between FM and Observation')
    ax = fig.add_subplot(111)
    ax.plot(ds_freq,y_obs,'r')
    ax.plot(ds_freq,y_FM,'b-',linewidth=2)
    #ax.set_ylim((-5,25))
    
    ax.set_xlabel('freq [Hz]')
    ax.set_ylabel('Tb [K]')
    
    ax.legend(('Observed','FM'))
    
    plt.show()
    pass

def retrieve_cycle(spectro_dataset, retrieval_param):
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
    print("Retrieval of Ozone and H20")

    cycle = retrieval_param["integration_cycle"]
    good_channels = spectro_dataset.good_channels[cycle].data == 1
    ds_freq = spectro_dataset.frequencies.values[good_channels]
    ds_y = spectro_dataset.Tb[cycle].values[good_channels]

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
    ac = arts.ArtsController(verbosity=2, agenda_verbosity=2)
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)

    # defining simulation grids
    if retrieval_param["binned_ch"]:
        f_grid = ds_freq
    else:
        f_grid = make_f_grid(retrieval_param)

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
    ac.set_spectroscopy_from_file(
        abs_lines_file = retrieval_param['line_file'],
        abs_species = ["O3","H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252", "O2-PWR93","N2-SelfContStandardType"],
        format = 'Arts',
        line_shape = ["VVH", 750e9],
        )
    
    #apriori_data_GROSOM.plot_apriori_cira86(retrieval_param)
    print('State of the atmosphere defined with:')
    if retrieval_param['atm'] == 'fascod':
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
    sensor = arts.SensorFFT(ds_freq, ds_df)
    ac.set_sensor(sensor)
    
    # doing the checks
    ac.checked_calc()
    
    # FM + noise --> to retrieve as test !
    y_FM = ac.y_calc()

    plot_FM_comparison(ds_freq,y_FM[0],ds_y)
    ac.set_y([ds_y])
    
    # Setup the retrieval
    z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
    z_top_ret = retrieval_param["z_top_ret_grid"]
    z_res_ret = retrieval_param["z_resolution_ret_grid"]
    z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
    p_grid_retrieval = z2p_simple(z_grid_retrieval)

    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])

    sx = covmat.covmat_diagonal_sparse(retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval))
    sx_water = covmat.covmat_diagonal_sparse(1e-5 * np.ones_like(p_grid_retrieval))
    
    ozone_ret = arts.AbsSpecies(
        species = 'O3',
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx,
        unit = 'vmr'
    )
    h2o_ret = arts.AbsSpecies(
        species = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx_water,
        unit = 'vmr'
    )

    #polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])
    
    # increase variance for spurious spectra by factor
    #factor = retrieval_param['increased_var_factor']
    
    #y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)
    
    #y_var[(level1b_ds.good_channels[cycle].values==0)] = factor*retrieval_param['unit_var_y']
    
    y_var = 15*spectro_dataset.stdTb[cycle].data[good_channels]
    #polyfit_ret = arts.Polyfit(
    #    poly_order=1, covmats=[np.array([[5]]), np.array([[1]])]
    #)

    ac.define_retrieval(retrieval_quantities=[ozone_ret, h2o_ret], y_vars=y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)
    
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

@arts_agenda
def inversion_iterate_agenda(ws):
    """Custom inversion iterate agenda to ignore bad partition functions."""
    ws.Ignore(ws.inversion_iteration_counter)

    # Map x to ARTS' variables
    ws.x2artsAtmAndSurf()
    ws.x2artsSensor()

    # To be safe, rerun some checks
    ws.atmfields_checkedCalc(negative_vmr_ok=1)
    ws.atmgeom_checkedCalc()

    # Calculate yf and Jacobian matching x
    ws.yCalc(y=ws.yf)

    # Add baseline term
    ws.VectorAddVector(ws.yf, ws.y, ws.y_baseline)

    # This method takes cares of some "fixes" that are needed to get the Jacobian
    # right for iterative solutions. No need to call this WSM for linear inversions.
    ws.jacobianAdjustAndTransform()

def retrieve_cycle_tropospheric_corrected(spectro_dataset, retrieval_param):
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
    ds_freq = spectro_dataset.frequencies.values[good_channels]
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
    ac = arts.ArtsController(verbosity=2, agenda_verbosity=2)
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)

    
    # defining simulation grids
    if retrieval_param["binned_ch"]:
        f_grid = ds_freq
    else:
        f_grid = make_f_grid(retrieval_param)

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
    print('State of the atmosphere defined with:')
    if retrieval_param['atm'] == 'fascod':
        atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param)
    elif retrieval_param['atm'] == 'fascod_2':
        ac.ws.AtmRawRead(basename = "planets/Earth/Fascod/{}/{}".format('midlatitude-winter','midlatitude-winter'))
        ac.ws.AtmFieldsCalc()        
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
    y_FM_noisy = y_FM+np.random.rand(len(y_FM[0]))
    
    plot_FM_comparison(ds_freq,y_FM[0],ds_y)

    ac.set_y([ds_y])
    
    # Setup the retrieval
    z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
    z_top_ret = retrieval_param["z_top_ret_grid"]
    z_res_ret = retrieval_param["z_resolution_ret_grid"]
    z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
    p_grid_retrieval = z2p_simple(z_grid_retrieval)

    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])
    
    #sx = covmat.covmat_diagonal_sparse(retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval))
    
    sx = covmat.covmat_1d_sparse(
        np.log10(p_grid_retrieval),
        retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval),
        0.3 * np.ones_like(p_grid_retrieval),
        fname="lin",
        cutoff=0.001,
    )

    sx = covmat.covmat_diagonal_sparse(retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval))
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
    #h2o_ret = arts.AbsSpecies('H2O', p_ret_grid, lat_ret_grid, lon_ret_grid, sx, unit='vmr')
    

    #polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])
    
    # increase variance for spurious spectra by factor
    #factor = retrieval_param['increased_var_factor']
    
    #y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)
    
    #y_var[(level1b_ds.good_channels[cycle].values==0)] = factor*retrieval_param['unit_var_y']
    
    #y_var = factor*utils.var_allan(ds_y) * np.ones_like(ds_y)
    y_var = 8*spectro_dataset.stdTb[cycle].data[good_channels]

    #polyfit_ret = arts.Polyfit(
    #    poly_order=1, covmats=[np.array([[5]]), np.array([[1]])]
    #)

    ac.define_retrieval(retrieval_quantities=[ozone_ret], y_vars=y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)
    
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

def test_retrieval(spectro_dataset, retrieval_param):
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
    print("Test of the retrieval setup")

    cycle = retrieval_param["integration_cycle"]
    good_channels = spectro_dataset.good_channels[cycle].data == 1
    ds_freq = spectro_dataset.frequencies.values[good_channels]
    ds_y = spectro_dataset.Tb[cycle].values[good_channels]
    #ds_y = spectro_dataset.Tb_corr[cycle].values[good_channels]

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
    ac = arts.ArtsController(verbosity=2, agenda_verbosity=2)
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)

    # defining simulation grids
    if retrieval_param["binned_ch"]:
        f_grid = ds_freq
    else:
        f_grid = make_f_grid(retrieval_param)

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
        abs_species=["O3", "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252", "O2-TRE05", "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252"],
        format='Arts',
        line_shape=["VVH", 750e9],
    )

    #apriori_data_GROSOM.plot_apriori_cira86(retrieval_param)
    print('State of the atmosphere defined with:')
    if retrieval_param['atm'] == 'fascod':
        atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param)
    elif retrieval_param['atm'] == 'fascod_2':
        ac.ws.AtmRawRead(basename = "planets/Earth/Fascod/{}/{}".format('midlatitude-winter','midlatitude-winter'))
        ac.ws.AtmFieldsCalc()
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
    y_FM_noisy = y_FM[0]+np.random.rand(len(y_FM[0]))

    plot_FM_comparison(ds_freq,y_FM_noisy,ds_y)
    ac.set_y([y_FM_noisy])
    
    # Setup the retrieval
    z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
    z_top_ret = retrieval_param["z_top_ret_grid"]
    z_res_ret = retrieval_param["z_resolution_ret_grid"]
    z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
    p_grid_retrieval = z2p_simple(z_grid_retrieval)

    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])

    sx = covmat.covmat_1d_sparse(
        np.log10(p_grid_retrieval),
        retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval),
        0.3 * np.ones_like(p_grid_retrieval),
        fname="lin",
        cutoff=0.001,
    )

    sx = covmat.covmat_diagonal_sparse(retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval))
    sx_water = covmat.covmat_diagonal_sparse(1e-13 * np.ones_like(p_grid_retrieval))
    
    ozone_ret = arts.AbsSpecies(
        species = 'O3',
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx,
        unit = 'vmr'
    )
    h2o_ret = arts.AbsSpecies(
        species = "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
        p_grid = p_grid_retrieval,
        lat_grid = lat_ret_grid,
        lon_grid = lon_ret_grid,
        covmat = sx_water,
        unit = 'vmr'
    )

    #polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])
    
    # increase variance for spurious spectra by factor
    #factor = retrieval_param['increased_var_factor']
    
    #y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)
    
    #y_var[(level1b_ds.good_channels[cycle].values==0)] = factor*retrieval_param['unit_var_y']
    
    y_var = 10*np.ones_like(y_FM_noisy)

    #polyfit_ret = arts.Polyfit(
    #    poly_order=1, covmats=[np.array([[5]]), np.array([[1]])]
    #)

    ac.define_retrieval(retrieval_quantities=[ozone_ret, h2o_ret], y_vars=y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)
    
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