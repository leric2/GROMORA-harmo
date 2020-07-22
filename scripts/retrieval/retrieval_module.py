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
from typhon.arts.xml import load

def make_f_grid(retrieval_param):
    '''
    create simulation frequency grid
    '''
    n_f = retrieval_param["number_of_freq_points"]# Number of points
    bw = 1.5e9  # Bandwidth
    x = np.linspace(-0.5, 0.5, n_f)
    f_grid = x ** 3 + x / 10
    f_grid = f_grid * bw / (max(f_grid) - min(f_grid))
    #f_grid = np.linspace(retrieval_param["f_min"]-10, retrieval_param["f_max"]+10, n_f)
    return f_grid

def plot_FM_comparison(ds_freq,f_grid,y_FM,y_obs):
    fig = plt.figure()
    
    fig.suptitle('Comparison between FM and Observation')
    ax = fig.add_subplot(111)
    ax.plot(ds_freq,y_obs,'r')
    ax.plot(ds_freq,y_FM[0],'b-',linewidth=2)
    #ax.set_ylim((-5,25))
    
    ax.set_xlabel('freq [Hz]')
    ax.set_ylabel('Tb [K]')
    
    ax.legend(('Observed','FM'))
    
    plt.show()
    pass

def retrieve_cycle(instrument, retrieval_param):
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
    level1b_ds = instrument.level1b_ds

    cycle = retrieval_param["integration_cycle"] 
    ds_freq = level1b_ds.frequencies.values
    ds_num_of_channel = len(ds_freq)
    ds_Tb = level1b_ds.Tb[cycle].values
    #ds_Tb_corr = level1b_ds.Tb_corr[cycle].values
    
    ds_bw = max(ds_freq) - min(ds_freq)
    
    ds_df = ds_bw/(ds_num_of_channel-1)
    
    retrieval_param["zenith_angle"]=level1b_ds.mean_sky_angle.values[cycle]
    retrieval_param["azimuth_angle"] = level1b_ds.azimuth_angle.values[cycle]

    retrieval_param["time"] = level1b_ds.time[cycle].values
    retrieval_param["lat"] = level1b_ds.lat[cycle].values
    retrieval_param["lon"] = level1b_ds.lon[cycle].values

    retrieval_param['time_start'] = level1b_ds.first_sky_time[cycle].values
    retrieval_param['time_stop'] = level1b_ds.last_sky_time[cycle].values
    
    retrieval_param["f_max"] = max(ds_freq)
    retrieval_param["f_min"] = min(ds_freq)
    retrieval_param["bandwidth"] = max(ds_freq) - min(ds_freq)
        
    # Iniializing ArtsController object
    ac = arts.ArtsController()
    
    # atmosphere + basic param
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)
    
    # defining simulation grids
    f_grid = make_f_grid(retrieval_param) + retrieval_param['obs_freq']
    z_bottom = retrieval_param["z_bottom_sim_grid"]
    z_top = retrieval_param["z_top_sim_grid"]
    z_res = retrieval_param["z_resolution_sim_grid"]
    z_grid = np.arange(z_bottom, z_top, z_res)
    p_grid = z2p_simple(z_grid)

    ac.set_grids(f_grid, p_grid)
    
    # altitude for the retrieval
    #ac.set_surface(level1b_ds.alt.values[cycle])
    ac.set_surface(1000)
    
    #spectroscopy
    ac.set_spectroscopy_from_file(
        abs_lines_file = retrieval_param['line_file'],
        abs_species = ["O3","H2O","H2O-PWR98", "O2-PWR98","N2-SelfContStandardType"],
        format = 'Arts',
        line_shape = ["VVH", 750e9],
        )
    
    t1 = pd.to_datetime(retrieval_param['time_start'])
    t2 = pd.to_datetime(retrieval_param['time_stop'])
    
    extra_time_ecmwf = 4

    #ecmwf_store = '/home/eric/Documents/PhD/GROSOM/ECMWF'
    ecmwf_store = retrieval_param['ecmwf_store_location'] 
    cira86_path = retrieval_param['cira86_path']

    atm = apriori_data_GROSOM.get_apriori_atmosphere_fascod_ecmwf_cira86(
        retrieval_param,
        ecmwf_store,
        cira86_path,
        t1,
        t2,
        extra_time_ecmwf
    )

    #ac.set_atmosphere_fascod('midlatitude-winter')
    #fascod_atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])
    ac.set_atmosphere(atm, vmr_zeropadding=True)

    # Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
    ac.apply_hse(100e2, 0.1)  # value taken from GROMOS retrieval
    
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
    ac.set_y([ds_Tb])
    
    # Defining our sensors
    sensor = arts.SensorFFT(ds_freq, ds_df)
    ac.set_sensor(sensor)
    
    # doing the checks
    ac.checked_calc()
    
    # FM
    y_FM = ac.y_calc()
    
    # Plotting a comparison to check which spectra is used
    plot_FM_comparison(ds_freq,f_grid,y_FM,ds_Tb)
    
    # Setup the retrieval
    z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
    z_top_ret = retrieval_param["z_top_ret_grid"]
    z_res_ret = retrieval_param["z_resolution_ret_grid"]
    z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
    p_grid_retrieval = z2p_simple(z_grid_retrieval)

    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])
    
    sx_O3 = covmat.covmat_diagonal_sparse(1e-6 * np.ones_like(p_grid_retrieval))
    sx_H2O = covmat.covmat_diagonal_sparse(1e-6 * np.ones_like(p_grid_retrieval))

    '''
    o3_std_value=0.7e-6
    o3_std_const = o3_std_value * np.ones_like(p_ret_grid)
    o3_apriori = atm.vmr_field('o3').data[:,0][:,0]
    o3_std_limit = p_interpolate(
        p_ret_grid, atm.vmr_field('o3').grids[0], o3_apriori, fill=0
    )

    o3_std = np.minimum(o3_std_limit, o3_std_const)
    std_correction_o3 = 1
    ozone_covmat = covmat.covmat_1d_sparse(
        np.log10(p_ret_grid),
        std_correction_o3 * o3_std,
        0.3 * np.ones_like(p_ret_grid),
        fname="lin",
        cutoff=0.001,
    )

    h2o_std_value=0.7e-6
    h2o_std_const = h2o_std_value * np.ones_like(p_ret_grid)
    h2o_apriori = atm.vmr_field('h2o').data[:,0][:,0]
    h2o_std_limit = p_interpolate(
        p_ret_grid, atm.vmr_field('h2o').grids[0], h2o_apriori, fill=0
    )

    h2o_std = np.minimum(h2o_std_limit, h2o_std_const)
    std_correction_h2o = 1
    h2o_covmat = covmat.covmat_1d_sparse(
        np.log10(p_ret_grid),
        std_correction_h2o * h2o_std,
        0.3 * np.ones_like(p_ret_grid),
        fname="lin",
        cutoff=0.001,
    )
    '''

    # The different things we want to retrieve
    ozone_ret = arts.AbsSpecies('O3', p_grid_retrieval, lat_ret_grid, lon_ret_grid, sx_O3, unit='vmr')
    h2o_ret = arts.AbsSpecies('H2O', p_grid_retrieval, lat_ret_grid, lon_ret_grid, sx_H2O, unit='vmr')
    
    #fshift_ret = arts.FreqShift(100e3, df=50e3)
    #polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])
    
    # increase variance for spurious spectra by factor
    
    factor = retrieval_param['increased_var_factor']
    
    
    y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)
    
    y_var[(level1b_ds.good_channels[cycle].values==0)] = factor*retrieval_param['unit_var_y']
    
    #y_var = utils.var_allan(ds_Tb) * np.ones_like(ds_Tb)

    ac.define_retrieval([ozone_ret, h2o_ret], y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)
    
    # Let a priori be off by 0.5 ppm (testing purpose)
    #vmr_offset = -0.5e-6
    #ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)
    
    # Run retrieval (parameter taken from MOPI)
    # SOMORA is using 'lm': Levenberg-Marquardt (LM) method
    ac.oem(
        method='lm',
        max_iter=10,
        stop_dx=0.01,
        lm_ga_settings=[100.0, 3.0, 5.0, 10.0, 1.0, 10.0],
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

def retrieve_cycle_tropospheric_corrected(instrument, retrieval_param):
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
    level1b_ds = instrument.level1b_ds
    #meteo_ds = instrument.meteo_ds

    cycle = retrieval_param["integration_cycle"] 
    ds_freq = level1b_ds.frequencies.values
    ds_num_of_channel = len(ds_freq)
    #ds_Tb = level1b_ds.Tb[cycle].values
    ds_Tb_corr = level1b_ds.Tb_corr[cycle].values
    
    ds_bw = max(ds_freq) - min(ds_freq)
    
    ds_df = ds_bw/(ds_num_of_channel-1)
    
    retrieval_param["zenith_angle"]=level1b_ds.mean_sky_angle.values[cycle]
    retrieval_param["azimuth_angle"] = level1b_ds.azimuth_angle.values[cycle]

    retrieval_param["time"] = level1b_ds.time[cycle].values
    retrieval_param["lat"] = level1b_ds.lat[cycle].values
    retrieval_param["lon"] = level1b_ds.lon[cycle].values

    retrieval_param['time_start'] = level1b_ds.first_sky_time[cycle].values
    retrieval_param['time_stop'] = level1b_ds.last_sky_time[cycle].values
    
    retrieval_param["f_max"] = max(ds_freq)
    retrieval_param["f_min"] = min(ds_freq)
    retrieval_param["bandwidth"] = max(ds_freq) - min(ds_freq)
        
    # Iniializing ArtsController object
    ac = arts.ArtsController()
    
    # atmosphere + basic param
    ac.setup(atmosphere_dim=1, iy_unit='RJBT', ppath_lmax=-1, stokes_dim=1)
    
    # defining simulation grids
    f_grid = make_f_grid(retrieval_param) + instrument.observation_frequency
    
    z_bottom = retrieval_param["z_bottom_sim_grid"]
    z_top = retrieval_param["z_top_sim_grid"]
    z_res = retrieval_param["z_resolution_sim_grid"]
    z_grid = np.arange(z_bottom, z_top, z_res)
    p_grid = z2p_simple(z_grid)

    ac.set_grids(f_grid, p_grid)
    
    # altitude for the retrieval
    #ac.set_surface(level1b_ds.alt.values[cycle])
    ac.set_surface(1500)
    
    # spectroscopy
    ac.set_spectroscopy_from_file(
        abs_lines_file = retrieval_param['line_file'],
        abs_species = ["O3","H2O-PWR98", "O2-PWR93","N2-SelfContStandardType"],
        format = 'Arts',
        line_shape = ("VVH", 750e9),
        )
    
    #apriori_data_GROSOM.plot_apriori_cira86(retrieval_param)

    #fascod_atm = apriori_data_GROSOM.get_apriori_fascod(retrieval_param)

    t1 = pd.to_datetime(retrieval_param['time_start'])
    t2 = pd.to_datetime(retrieval_param['time_stop'])
    
    extra_time_ecmwf = 4

    #ecmwf_store = '/home/eric/Documents/PhD/GROSOM/ECMWF'
    ecmwf_store = retrieval_param['ecmwf_store_location'] 
    cira86_path = retrieval_param['cira86_path']

    atm = apriori_data_GROSOM.get_apriori_atmosphere_fascod_ecmwf_cira86(
        retrieval_param,
        ecmwf_store,
        cira86_path,
        t1,
        t2,
        extra_time_ecmwf
    )

    #ac.set_atmosphere_fascod('midlatitude-winter')
    #fascod_atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])
    ac.set_atmosphere(atm, vmr_zeropadding=True)

    # Apply HSE, only allowed after setup of the atmosphere, when z is already on simulation grid
    ac.apply_hse(100e2, 0.1)  # value taken from GROMOS retrieval

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
    ac.set_y([ds_Tb_corr])
    
    # Defining our sensors
    sensor = arts.SensorFFT(ds_freq, ds_df)
    ac.set_sensor(sensor)
    
    # doing the checks
    ac.checked_calc()
    
    # FM
    y = ac.y_calc()
    
    plot_FM_comparison(ds_freq,f_grid,y,ds_Tb_corr)

    # Setup the retrieval
    z_bottom_ret = retrieval_param["z_bottom_ret_grid"]
    z_top_ret = retrieval_param["z_top_ret_grid"]
    z_res_ret = retrieval_param["z_resolution_ret_grid"]
    z_grid_retrieval = np.arange(z_bottom_ret, z_top_ret, z_res_ret)
    p_grid_retrieval = z2p_simple(z_grid_retrieval)

    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])
    
    sx = covmat.covmat_diagonal_sparse(retrieval_param["apriori_O3_cov"] * np.ones_like(p_grid_retrieval))
    
    # The different things we want to retrieve
    ozone_ret = arts.AbsSpecies('O3', p_grid_retrieval, lat_ret_grid, lon_ret_grid, sx, unit='vmr')
    #h2o_ret = arts.AbsSpecies('H2O', p_ret_grid, lat_ret_grid, lon_ret_grid, sx, unit='vmr')
    
    #fshift_ret = arts.FreqShift(100e3, df=50e3)
    #polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])
    
    # increase variance for spurious spectra by factor
    factor = retrieval_param['increased_var_factor']
    
    y_var = retrieval_param['unit_var_y'] * np.ones_like(ds_freq)
    
    y_var[(level1b_ds.good_channels[cycle].values==0)] = factor*retrieval_param['unit_var_y']
    
    #y_var = utils.var_allan(ds_Tb) * np.ones_like(ds_Tb)

    ac.define_retrieval([ozone_ret], y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)
    
    # Let a priori be off by 0.5 ppm (testing purpose)
    #vmr_offset = -0.5e-6
    #ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)
    
    # Run retrieval (parameter taken from MOPI)
    # SOMORA is using 'lm': Levenberg-Marquardt (LM) method
    ac.oem(
        method='gn',
        max_iter=5,
        stop_dx=0.1,
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