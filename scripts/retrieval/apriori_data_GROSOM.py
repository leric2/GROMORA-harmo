#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 

@author: eric

Collection of functions for dealing with GROSOM a-priori data

Using the Atmosphere Class to define a 

Including : 
    * a-priori data
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt

from retrievals import arts
from retrievals.data.ecmwf import ECMWFLocationFileStore
from retrievals.data.ecmwf import levels
from retrievals.data import p_interpolate

from typhon.arts.xml import load

class APrioriDataGROSOM(arts.Atmosphere):
    '''
    TO DO ?
    '''
    def __init__(self, apriori_type):
        self.type = apriori_type

        pass

def extract_ecmwf_ds(ECMWF_store_path, t1, t2):
    '''
    Building the ecmwf store for atmospheric state
    '''
    ecmwf_store = ECMWFLocationFileStore('/home/eric/Documents/PhD/GROSOM/ECMWF', 'ecmwf_oper_v2BERN_%Y%m%d.nc')
    ecmwf_ds = (
        ecmwf_store.select_time(t1, t2 )
        .mean(dim='time')
        .swap_dims({"level":"pressure"})
    )
    return ecmwf_ds

def plot_apriori_cira86(retrieval_param):
    lat = retrieval_param['lat']
    month = pd.to_datetime(retrieval_param['time']).month

    cira86 = read_cira86_monthly(retrieval_param['cira86_path'], month, lat)

    plot_cira86_profile(cira86)
    pass 

def get_apriori_fascod(retrieval_param):
    fascod_atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])
    return fascod_atm

def get_apriori_atmosphere(retrieval_param):
    '''
    Defining a-priori atmosphere from ECMWF operationnal dataset and Fascod clim

    First, we initiate an Atmosphere from Fascod climatology (not needed) before
    updating the required fields from ECMWF operationnal dataset.

    NOT WORKING YET
    '''
    
    ECMWF_store_path = retrieval_param['ecmwf_store_path']

    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, t1 = '2019-04-16 03:30', t2 = '2019-04-16 16:45')

    atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])

    # Temperature
    atm.set_t_field(ds_ecmwf['pressure'].values, ds_ecmwf['temperature'].values)

    # z field
    atm.set_z_field(ds_ecmwf["pressure"].values, ds_ecmwf["pressure"].values)

    # Ozone
    atm.set_vmr_field(
        "O3", ds_ecmwf["pressure"].values, ds_ecmwf['ozone_mass_mixing_ratio'].values
        )
    
    return atm

def get_apriori_atmosphere_ecmwf_cira86(retrieval_param):
    '''
    Defining a-priori atmosphere from ECMWF operationnal dataset and CIRA86 clim

    First, we initiate an Atmosphere from Fascod climatology (not needed) before
    updating the required fields from ECMWF operationnal dataset.

    TO DO: add the CIRA86 cimatology, otherwise the ecmwf data do not go high enough.
    '''
    
    # Read EACMWF data (oper for now)
    ECMWF_store_path = retrieval_param['ecmwf_store_path']
    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, t1 = '2019-04-16 03:30', t2 = '2019-04-16 16:45')
    
    # Reading CIRA86
    lat = retrieval_param['lat']
    month = pd.to_datetime(retrieval_param['time']).month

    cira86 = read_cira86_monthly(retrieval_param['cira86_path'], month, lat)

    #plot_ecmwf_cira86_profile(ds_ecmwf, cira86)

    # Merging ecmwf and CIRA86
    ds_ptz = merge_ecmwf_cira86(ds_ecmwf, cira86)

    
    plot_apriori_ptz(ds_ptz)

    '''
    At this point we have 2 choice:
    1. Create atm using the fascod and overwrite our variable from ds_ptz

    2. Create atm using "from_dataset". 
    In this case:
    The dataset must have the coordinates `pressure`, `lat`, `lon` and can have
    variables `t`, `z` and absorbing species like `o3` (all lowercase).
    --> It means we have to expand_dims at that point to use this method.
    --> It means also that we need apriori values for each of the species in 
    abs_species ......

    '''

    #ds_ptz = ds_ptz.expand_dims(dim=['lat','lon'])
    # atm = arts.Atmosphere.from_dataset(ds_ptz)
    atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])

    # Temperature
    atm.set_t_field(ds_ptz['p'].values, ds_ptz['t'].values)

    # z field
    atm.set_z_field(ds_ptz["p"].values, ds_ptz["z"].values)
    
    # Ozone
    atm.set_vmr_field(
        "O3", ds_ecmwf["pressure"].values, ds_ecmwf['ozone_mass_mixing_ratio'].values
        )

    # Water vapor
    atm.set_vmr_field(
        "H2O", ds_ecmwf["pressure"].values, ds_ecmwf['specific_humidity'].values
        )
    
    return atm

def merge_ecmwf_cira86(ds_ecmwf, cira86):
    '''
    Merging profile from ECMWF oper and CIRA86 monthly climatology.

    Start with a very simple scheme, we just take CIRA86 values when we reach the top
    of the ECMWF ones.

    Parameters
        ----------
        ds_ecmwf : xarray.Dataset
            
        cira86 : xarray.Dataset
    
    Returns
        -------
        TYPE
            DESCRIPTION.
    '''

    upper_p_grid = cira86.Pressure.data[cira86.Pressure < np.min(ds_ecmwf.pressure.data)]
    #upper_t = cira86.temperature.data[cira86.Pressure < np.min(ds_ecmwf.pressure.data)]
    
    upper_cira86_ds = cira86.sel(Pressure = upper_p_grid)

    p_grid = np.hstack((ds_ecmwf.pressure.data, upper_p_grid))

    #ecmwf_t_i = p_interpolate(
    #    p_grid, ds_ecmwf.pressure.data, ds_ecmwf.temperature.data, fill=np.nan
    #)

    temperature = np.hstack((ds_ecmwf.temperature.data, upper_cira86_ds.temperature.data))

    # interpolate CIRA86 altitude on p_grid
    cira86_alt_i = p_interpolate(
        p_grid, cira86.Pressure.data, cira86.altitude.data, fill = np.nan
    )

    ds_merged = xr.Dataset({'t': ('p', temperature),
                            'z': ('p', cira86_alt_i)},
                            coords = {
                                'p' : ('p', p_grid),
                            }
    )

    return ds_merged

def read_cira86_monthly(cira86_datapath, month = None, latitude = None):
    '''
    Extract profile from the CIRA86 monthly climatology.

    For now, works with the ones saved in arts-data saved as .xml files:
    2 files per month in the form:
    cira86_month*.t.xml -> temperature data, zonally averaged
    cira86_month*.z.xml -> geometric altitude


    Parameters
        ----------
        month : INT
            
        latitude : INT
    
    Returns
        -------
        TYPE
            DESCRIPTION.
    '''
    
    cira86_filename_t = 'cira86_month{}.t.xml'.format(month)
    cira86_filename_z = 'cira86_month{}.z.xml'.format(month)

    cira86_t = load(os.path.join(cira86_datapath, cira86_filename_t)).to_xarray()
    cira86_z = load(os.path.join(cira86_datapath, cira86_filename_z)).to_xarray()

    # merge the 2 ds
    cira86 = xr.merge([cira86_z.rename('altitude'),cira86_t.rename('temperature')], compat='equals', join = 'inner')

    # Select the right latitude
    cira86_lat = np.round(latitude, -1)

    if cira86_lat not in cira86_t.Latitude:
        raise ValueError('No valid latitude')

    return cira86.sel(Latitude = cira86_lat, Longitude = 0)

def plot_ecmwf_cira86_profile(ds_ecmwf, cira86):
    fig, axs = plt.subplots(1,2, sharex= True, sharey = True)
    axs[0].plot(cira86.temperature, cira86.Pressure, label='Temperature')
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('T [K]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(ds_ecmwf.temperature, ds_ecmwf.pressure, label='Temperature')
    #axs[1].invert_yaxis()
    axs[1].set_yscale('log')
    axs[1].set_xlabel('T [K]')
    axs[1].set_ylabel('$P$ [Pa]')

    #fig.suptitle('CIRA86 Temperature profile')
    fig.tight_layout()

    fig.show()
    pass

def plot_apriori_ptz(ds_ptz):
    fig, axs = plt.subplots(1,2, sharex= True)
    axs[0].plot(ds_ptz.t, ds_ptz.p)
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('T [K]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(ds_ptz.t, ds_ptz.z/1e3)
    #axs[1].invert_yaxis()
    #axs[1].set_yscale('log')
    axs[1].set_xlabel('T [K]')
    axs[1].set_ylabel('$Z$ [km]')

    #fig.suptitle('CIRA86 Temperature profile')
    fig.tight_layout()

    fig.show()
    pass

def plot_cira86_profile(cira86):
    fig, axs = plt.subplots(1,2)
    axs[0].plot(cira86.temperature, cira86.Pressure, label='Temperature')
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('T [K]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(cira86.temperature, cira86.altitude, label='Temperature')
    axs[1].set_xlabel('T [K]')
    axs[1].set_ylabel('$Z$ [m]')

    fig.suptitle('CIRA86 Temperature profile')
    fig.tight_layout()

    fig.show()
    pass