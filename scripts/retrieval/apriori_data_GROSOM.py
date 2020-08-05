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

ARTS_DATA_PATH = os.environ.get("ARTS_DATA_PATH", None)

class APrioriDataGROSOM(arts.Atmosphere):
    '''
    TO DO ?
    '''
    def __init__(self, apriori_type):
        self.type = apriori_type

        pass

def extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix, t1, t2):
    '''i
    Building the ecmwf store for atmospheric state
    '''
    ecmwf_store = ECMWFLocationFileStore(ECMWF_store_path, ecmwf_prefix)
    ecmwf_ds = (
        ecmwf_store.select_time(t1, t2, combine='by_coords')
        .mean(dim='time')
        .swap_dims({"level":"pressure"})
    )
    return ecmwf_ds

def read_o3_apriori_ecmwf_mls_gromosOG(filename):
    '''
    read the apriori o3 used in gromos retrieval
    
    just for fun, no idea where that file come from...
    '''
    ds = pd.read_csv(
        filename,
        sep=' ',
        skiprows=6,
        header=None,
        usecols=[0,1]
        )

    o3_apriori_ds = xr.Dataset({'o3': ('p', ds.iloc[:,1])},
                            coords = {
                                'p' : ('p', ds.iloc[:,0]),
                            }
    )
    return o3_apriori_ds

def read_o3_apriori_OG_SOMORA(filename, m):
    '''
    read the apriori o3 used in somora retrieval
    
    just for fun, no idea where that file come from...

    is a mix of MLS and RS ?
    '''
    ds = pd.read_csv(
        filename,
        sep=',',
        skiprows=0,
        header=None,
        #usecols=[0,1]
        )

    o3_apriori_ds = xr.Dataset({'o3': (['altitude','month'], ds.iloc[:,1:])},
                            coords = {
                                'altitude' : ('altitude', ds.iloc[:,0]),
                                'month' : ('month', np.arange(1,13,1))
                            }
    )
    return o3_apriori_ds.sel(month = m)

def plot_apriori_cira86(retrieval_param):
    lat = retrieval_param['lat']
    month = pd.to_datetime(retrieval_param['time']).month

    cira86 = read_cira86_monthly(retrieval_param['cira86_path'], month, lat)

    plot_cira86_profile(cira86)
    pass 

def get_apriori_fascod(retrieval_param):
    '''
    All in one.
    '''
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

    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix, t1 = '2019-04-16 03:30', t2 = '2019-04-16 16:45')

    atm = arts.Atmosphere.from_arts_xml(retrieval_param['prefix_atm'])

    # Temperature
    atm.set_t_field(ds_ecmwf['pressure'].values, ds_ecmwf['temperature'].values)

    # z field
    atm.set_z_field(ds_ecmwf["pressure"].values, ds_ecmwf["pressure"].values)

    # Ozone
    #atm.set_vmr_field(
    #    "O3", ds_ecmwf["pressure"].values, ds_ecmwf['ozone_mass_mixing_ratio'].values
    #    )
    
    return atm

def compare_ecmwf_oper_dataset():
    # for comparing between ECMWF oper version or between GRIB and netCDF file
    ECMWF_store_path = '/home/esauvageat/Documents/GROSOM/Analysis/ECMWF'
    location = 'BERN'

    t1 = pd.to_datetime('2020-01-03 03:00:00')
    t2 = pd.to_datetime('2020-01-03 09:00:00')

    ecmwf_prefix = f'ecmwf_oper_v{2}_{location}_%Y%m%d.nc'
    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix, t1, t2)
    
    ecmwf_prefix2 = f'ecmwf_oper_v{2}_from_grib{location}_%Y%m%d.nc'
    ds_ecmwf_2 = extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix2, t1, t2)
    
    plot_ecmwf_comparison(ds_ecmwf, ds_ecmwf_2)

def get_apriori_atmosphere_fascod_ecmwf_cira86(retrieval_param, ecmwf_store, cira86_path, t1, t2, extra_time_ecmwf):
    '''
    Defining a-priori atmosphere from ECMWF operationnal dataset and CIRA86 clim

    First, we initiate an Atmosphere from Fascod climatology before
    updating the fields from ECMWF operationnal dataset.

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
    summer_months = [4, 5, 6, 7, 8, 9]

    lat = retrieval_param['lat']
    month = pd.to_datetime(retrieval_param['time']).month

    if month in summer_months:
        fascod_clim = 'midlatitude-summer'
    else:
        fascod_clim = 'midlatitude-winter'

    clim_prefix = os.path.join(ARTS_DATA_PATH, "planets/Earth/Fascod/{0}/{0}.".format(fascod_clim))
    
    # atm = arts.Atmosphere.from_dataset(ds_ptz)
    atm = arts.Atmosphere.from_arts_xml(clim_prefix)

    ecmwf_time1 = t1 - pd.Timedelta(extra_time_ecmwf/2, "h")
    ecmwf_time2 = t2 + pd.Timedelta(extra_time_ecmwf/2, "h")

    # Read EACMWF data (oper for now)
    ECMWF_store_path = ecmwf_store
    location = 'BERN'
    ecmwf_prefix = f'ecmwf_oper_v{2}_{location}_%Y%m%d.nc'
    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix, ecmwf_time1, ecmwf_time2)

    if sum(np.isnan(ds_ecmwf.pressure.values)) > 1:
        raise ValueError('no ECMWF data')

    # Reading CIRA86
    cira86 = read_cira86_monthly(cira86_path, month, lat)

    plot_ecmwf_cira86_profile(ds_ecmwf, cira86)
    

    o3_apriori = read_o3_apriori_ecmwf_mls_gromosOG(retrieval_param['apriori_ozone_climatology_GROMOS'])

    #filename = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
    
    o3_apriori_SOMORA = read_o3_apriori_OG_SOMORA(retrieval_param['apriori_ozone_climatology_SOMORA'], month)

    # Merging ecmwf and CIRA86
    ds_ptz = merge_ecmwf_cira86(ds_ecmwf, cira86)
    plot_apriori_ptz(ds_ptz)
    
    # Temperature
    atm.set_t_field(ds_ptz['p'].values, ds_ptz['t'].values)

    # z field
    atm.set_z_field(ds_ptz["p"].values, ds_ptz["z"].values)
    
    # DO NOT ADD O3 from ECMWF --> no value over 2 Pa...
    # Ozone
    compare_o3_apriori_OG(o3_apriori,o3_apriori_SOMORA)
    atm.set_vmr_field(
        "O3", o3_apriori["p"].values, o3_apriori['o3'].values
        )
    
    # Water vapor
    # merging Fascod with ECMWF
    #h2o_apriori = merge_ecmwf_Fascod_atm(ds_ecmwf, atm)
    #plot_apriori_ptz(h2o_apriori)
    
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

def merge_ecmwf_Fascod_atm(ds_ecmwf, fascod_atm):
    '''
    Merging profile from ECMWF oper and Fascod monthly climatology.

    Start with a very simple scheme, we just take Fascod values when we reach the top
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
    fascod = fascod_atm.vmr_field('h2o').to_xarray()
    fascod_altitude = fascod_atm.z_field.to_xarray()
    fascod = fascod.sel(Latitude=0, Longitude=0)
    fascod_altitude = fascod_altitude.sel(Latitude=0, Longitude=0)

    upper_p_grid = fascod.Pressure.data[fascod.Pressure < np.min(ds_ecmwf.pressure.data)]
    #upper_t = cira86.temperature.data[cira86.Pressure < np.min(ds_ecmwf.pressure.data)]
    
    upper_Fascod_ds = fascod.sel(Pressure = upper_p_grid)

    p_grid = np.hstack((ds_ecmwf.pressure.data, upper_p_grid))

    #ecmwf_t_i = p_interpolate(
    #    p_grid, ds_ecmwf.pressure.data, ds_ecmwf.temperature.data, fill=np.nan
    #)

    specific_humidity = np.hstack((ds_ecmwf.specific_humidity.data, upper_Fascod_ds.data))

    # interpolate Fascod altitude on p_grid
    fascod_alt_i = p_interpolate(
        p_grid, fascod.Pressure.data, fascod_altitude.data, fill = np.nan
    )

    ds_merged = xr.Dataset({'t': ('p', specific_humidity),
                            'z': ('p', fascod_alt_i)},
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

    fig.suptitle('PTZ profile from ECMWF and CIRA86')

    fig.show()

def plot_ecmwf_comparison(ds_ecmwf, ds_ecmwf_2):

    print(f'Geopotential :{ds_ecmwf.geopotential.values} vs {ds_ecmwf_2.geopotential.values}')
    print()
    print(f'lnsp :{ds_ecmwf.logarithm_of_surface_pressure.values} vs {ds_ecmwf_2.logarithm_of_surface_pressure.values}')
    print()

    #plt.rcParams['axes.grid'] = True

    fig, axs = plt.subplots(1,4, sharex= False, sharey = True)
    
    axs[0].plot(ds_ecmwf.temperature, ds_ecmwf.pressure, 'k')
    axs[0].plot(ds_ecmwf_2.temperature, ds_ecmwf_2.pressure, 'r--')
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('T [K]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(ds_ecmwf.ozone_mass_mixing_ratio*1e6, ds_ecmwf.pressure, 'k')
    axs[1].plot(ds_ecmwf_2.ozone_mixing_ratio*1e6, ds_ecmwf_2.pressure, 'r--')
    #axs[1].invert_yaxis()
    axs[1].set_yscale('log')
    axs[1].set_xlabel('O3 [ppm]')
    axs[1].set_ylabel('$P$ [Pa]')

    axs[2].plot(ds_ecmwf.specific_humidity, ds_ecmwf.pressure, 'k')
    axs[2].plot(ds_ecmwf_2.specific_humidity, ds_ecmwf_2.pressure, 'r--')
    #axs[1].invert_yaxis()
    axs[2].set_yscale('log')
    axs[2].set_xlabel('specific humidity')
    axs[2].set_ylabel('$P$ [Pa]')

    axs[3].plot(ds_ecmwf.vorticity_relative, ds_ecmwf.pressure, 'k')
    axs[3].plot(ds_ecmwf_2.relative_vorticity, ds_ecmwf_2.pressure, 'r--')
    #axs[1].invert_yaxis()
    axs[3].set_yscale('log')
    axs[3].set_xlabel('relative vorticity')
    axs[3].set_ylabel('$P$ [Pa]')
    
    fig.suptitle('Comparison between ECMWF datasets')

    fig.show()

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

    fig.suptitle('Apriori ptz profile')

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

    fig.show()
    pass

def plot_o3_apriori(o3_apriori):
    fig, ax = plt.subplots()
    ax.plot(o3_apriori.o3*1e6, o3_apriori.p)
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_xlabel('$O_3$ [PPM]')
    ax.set_ylabel('$P$ [Pa]')

    fig.suptitle('$O_3$ apriori from GROMOS OG')

    fig.show()
    pass

def compare_o3_apriori_OG(o3_apriori, o3_apriori_SOMORA):
    fig, axs = plt.subplots(1,2, sharex=True)
    axs[0].plot(o3_apriori.o3*1e6, o3_apriori.p)
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('$O_3$ [PPM]')
    axs[0].set_ylabel('$P$ [Pa]')

    axs[1].plot(o3_apriori_SOMORA.o3, o3_apriori_SOMORA.altitude)
    axs[1].set_xlabel('$O_3$ [PPM]')
    axs[1].set_ylabel('$Z$ [km]')

    fig.suptitle('$O_3$ apriori from OG retrievals')

    fig.show()