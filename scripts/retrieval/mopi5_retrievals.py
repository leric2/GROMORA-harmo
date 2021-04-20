#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:58:02 2020

@author: eric

Module with function for the retrieval of Ozone from MOPI5 using pyretrievals

"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

import retrieval_module
import apriori_data_GROSOM
#import GROSOM_library

from retrievals import arts
from retrievals import covmat
from retrievals import utils

from retrievals.arts.atmosphere import p2z_simple, z2p_simple
from retrievals.data import p_interpolate

# from typhon.arts.workspace import arts_agenda
from pyarts.workspace import arts_agenda

from utils_GROSOM import var_allan

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

def retrieve_cycle_tropospheric_corrected_mopi5(spectro_dataset, retrieval_param):
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
    noise_level = spectro_dataset.noise_level[cycle].values

    retrieval_param["zenith_angle"] = retrieval_param['ref_elevation_angle'] - spectro_dataset.mean_sky_elevation_angle.values[cycle]
    retrieval_param["azimuth_angle"] = spectro_dataset.azimuth_angle.values[cycle]
    
    retrieval_param["lat"] = spectro_dataset.lat[cycle].values
    retrieval_param["lon"] = spectro_dataset.lon[cycle].values
    try:
        retrieval_param["time"] = spectro_dataset.time[cycle].values
        retrieval_param['time_start'] = spectro_dataset.first_sky_time[cycle].values
        retrieval_param['time_stop'] = spectro_dataset.last_sky_time[cycle].values
    except:
        retrieval_param["time"] = 1
        retrieval_param['time_start'] = 1
        retrieval_param['time_stop'] = 2
    retrieval_param["f_max"] = max(ds_freq)
    retrieval_param["f_min"] = min(ds_freq)
    retrieval_param["bandwidth"] = 1e9

    ds_num_of_channel = len(ds_freq)
    #ds_Tb = Tb[cycle].values
    
    ds_bw = max(ds_freq) - min(ds_freq)
    
    ds_df = ds_bw/(ds_num_of_channel-1)
        
    # Iniializing ArtsController object
    ac = arts.ArtsController(verbosity=0, agenda_verbosity=0)
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
        time = cycle
        )

    ac.set_observations([obs])
    
    # Defining our sensors
    sensor = arts.SensorFFT(ds_freq, ds_df)
    ac.set_sensor(sensor)
    
    # doing the checks
    ac.checked_calc()
    
    # FM + noise --> to retrieve as test !
    y_FM = ac.y_calc()
    
    retrieval_module.plot_FM_comparison(ds_freq,y_FM[0],ds_y)

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
    
    #y_var[(level1b_ds.good_channels[cycle].values==0)] = factor*retrieval_param['unit_var_y']
    
    #y_var = 10*spectro_dataset.noise_level[cycle].data * np.ones_like(ds_y)
    y_var = retrieval_param['increased_var_factor']*np.square(spectro_dataset.stdTb[cycle].data[good_channels])
    fig, ax = plt.subplots(1,1,sharex=True)
    y_var = 100*var_allan(ds_y) * np.ones(len(ds_y))
    ax.plot(y_var, label='allan_var')
    y_var = retrieval_param['increased_var_factor']*np.square(noise_level) * np.ones(len(ds_y))
    ax.plot(y_var, label='noise_level')
    ax.legend()

    polyfit_ret = arts.Polyfit(
        poly_order=1, covmats=[np.array([[5]]), np.array([[1]])]
    )

    #ac.define_retrieval(retrieval_quantities=[ozone_ret], y_vars=y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)
    ac.define_retrieval_compatiblepyarts(retrieval_quantities=[ozone_ret], y_vars=y_var)
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

if __name__=="__main__":
    pass