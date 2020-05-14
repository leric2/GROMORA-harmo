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

from retrievals import arts
from retrievals import covmat
from retrievals import utils
from typhon.arts.workspace import arts_agenda


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

def plot_FM_comparison(ds_freq,f_grid,y,ds_Tb_corr):
    fig= plt.figure()

    ax = fig.add_subplot(111)
    ax.plot(ds_freq,ds_Tb_corr,'r')
    ax.plot(ds_freq,y[0],'b-',linewidth=2)
    ax.set_ylim((-5,25))
    pass

def plot(ds, ac, retrieval_param, title=""):
    '''
    Plotting function directly taken from Jonas ;)
    OG can be found in retrieval.py in MOPI retrievals
    
    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    ac : TYPE
        DESCRIPTION.
    title : TYPE, optional
        DESCRIPTION. The default is "".

    Returns
    -------
    figures : TYPE
        DESCRIPTION.

    '''
    ozone_ret, = ac.retrieval_quantities

    f_backend = ds.frequencies.values
    y = ds.Tb_trop_corr[retrieval_param['integration_cycle']].values
    yf = ac.yf[0]
    r = y - yf
    r_smooth = np.convolve(r, np.ones((128,)) / 128, mode="same")

    figures = list()

    fig, axs = plt.subplots(2, sharex=True)
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, y, label="observed")
    axs[0].plot((f_backend - retrieval_param['obs_freq']) / 1e6, yf, label="fitted")
    axs[0].set_ylim(-15, 50)
    axs[0].legend()
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r, label="residuals")
    axs[1].plot((f_backend - retrieval_param['obs_freq']) / 1e6, r_smooth, label="residuals smooth")
    axs[1].set_ylim(-15, 15)
    axs[1].legend()
    axs[1].set_xlabel("f - {:.3f} GHz [MHz]".format(retrieval_param['obs_freq'] / 1e9))

    for ax in axs:
        ax.set_ylabel("$T_B$ [K]")
        ax.set_xlim([-500, 700])
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)

    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].plot(
        ozone_ret.x * 1e6, ozone_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].plot(ozone_ret.xa * 1e6, ozone_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlabel("Ozone VMR [ppm]")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(ozone_ret.mr, ozone_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")

    axs[1][0].plot(ozone_ret.eo * 1e6, ozone_ret.z_grid / 1e3)
    axs[1][0].set_xlabel("$e_o$ [ppm]")
    axs[1][0].set_ylabel("Altitude [km]")

    for avk in ozone_ret.avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1][1].plot(avk, ozone_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")
    fig.suptitle(title + " Ozone")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)
    
    '''
    fig, axs = plt.subplots(2, 2, sharey=True)
    axs[0][0].plot(
        h2o_ret.x * 1e6, h2o_ret.z_grid / 1e3, label="retrieved", marker="x"
    )
    axs[0][0].plot(h2o_ret.xa * 1e6, h2o_ret.z_grid / 1e3, label="apriori")
    axs[0][0].set_xlabel("Water VMR [ppm]")
    axs[0][0].set_ylabel("Altitude [km]")
    axs[0][0].legend()

    axs[0][1].plot(h2o_ret.mr, h2o_ret.z_grid / 1e3)
    axs[0][1].set_xlabel("Measurement response")

    axs[1][0].plot(h2o_ret.eo * 1e6, h2o_ret.z_grid / 1e3)
    axs[1][0].set_xlabel("$e_o$ [ppm]")
    axs[1][0].set_ylabel("Altitude [km]")

    for avk in h2o_ret.avkm:
        if 0.8 <= np.sum(avk) <= 1.2:
            axs[1][1].plot(avk, h2o_ret.z_grid / 1e3)
    axs[1][1].set_xlabel("AVKM")
    fig.suptitle(title + " Water (v{})".format(1))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    figures.append(fig)
    '''
    return figures

def retrieve_cycle(level1b_dataset, meteo_ds, retrieval_param):
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
    ds_freq = level1b_dataset.frequencies.values
    ds_num_of_channel = len(ds_freq)
    ds_Tb = level1b_dataset.Tb[cycle].values
    ds_Tb_corr = level1b_dataset.Tb_trop_corr[cycle].values
    
    ds_bw = max(ds_freq) - min(ds_freq)
    
    ds_df = ds_bw/(ds_num_of_channel-1)
    
    retrieval_param["zenith_angle"]=level1b_dataset.meanAngleAntenna.values[cycle]
    
    retrieval_param["time"] = level1b_dataset.time[cycle].values
    retrieval_param["lat"] = level1b_dataset.lat[cycle].values
    retrieval_param["lon"] = level1b_dataset.lon[cycle].values
    
    retrieval_param["f_max"] = max(ds_freq)
    retrieval_param["f_min"] = min(ds_freq)
    retrieval_param["bandwidth"] = max(ds_freq) - min(ds_freq)
        
    # Iniializing ArtsController object
    ac = arts.ArtsController()
    
    # atmosphere + basic param
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
    #fascod_atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])
    #ac.set_atmosphere(fascod_atm)
    
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
    y = ac.y_calc()
    
    plot_FM_comparison(ds_freq,f_grid,y,ds_Tb_corr)
    
    # Setup the retrieval
    p_ret_grid = np.logspace(5, -1, 81)
    lat_ret_grid = np.array([retrieval_param["lat"]])
    lon_ret_grid = np.array([retrieval_param["lon"]])
    
    sx = covmat.covmat_diagonal_sparse(1e-6 * np.ones_like(p_ret_grid))
    
    # The different things we want to retrieve
    ozone_ret = arts.AbsSpecies('O3', p_ret_grid, lat_ret_grid, lon_ret_grid, sx, unit='vmr')
    #h2o_ret = arts.AbsSpecies('H2O', p_ret_grid, lat_ret_grid, lon_ret_grid, sx, unit='vmr')
    
    #fshift_ret = arts.FreqShift(100e3, df=50e3)
    #polyfit_ret = arts.Polyfit(poly_order=1, covmats=[np.array([[0.5]]), np.array([[1.4]])])

    y_var = 2e-2 * np.ones_like(ds_freq)
    
    #y_var = utils.var_allan(ds_Tb) * np.ones_like(ds_Tb)

    ac.define_retrieval([ozone_ret], y_var)
    #ac.define_retrieval([ozone_ret, fshift_ret, polyfit_ret], y_var)

    # Let a priori be off by 0.5 ppm
    #vmr_offset = 0.5e-6
    #ac.ws.Tensor4AddScalar(ac.ws.vmr_field, ac.ws.vmr_field, vmr_offset)
    
    # Run retrieval
    ac.oem(
        method='gn',
        max_iter=10,
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