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
import sys
import numpy as np
import xarray as xr
import pandas as pd
import math
import netCDF4
import matplotlib.pyplot as plt
import datetime as dt
from pytz import timezone


from retrievals import arts
from retrievals.data.ecmwf import ECMWFLocationFileStore
from retrievals.data.ecmwf import levels
from retrievals.data import interpolate
from retrievals.data import p_interpolate

from GROMORA_time import pysolar_sza, get_LST_from_GROMORA

#from typhon.arts.xml import load
from pyarts.xml import load

ARTS_DATA_PATH = os.environ.get("ARTS_DATA_PATH", None)

summer_months = [4, 5, 6, 7, 8, 9]

month_start = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
month_stop = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

class APrioriDataGROSOM(arts.Atmosphere):
    '''
    TO DO ?
    '''
    def __init__(self, apriori_type):
        self.type = apriori_type

        pass

def extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix, t1, t2):
    '''
    Building the ecmwf store for atmospheric state
    '''
    if t1 > dt.datetime(2013,6,24):
        ecmwf_levels=137
    elif t1 < dt.datetime(2013,6,25):
        ecmwf_levels=91
    else:
        ecmwf_levels=np.nan
        
    ecmwf_store = ECMWFLocationFileStore(ECMWF_store_path, ecmwf_prefix)
    ds_ecmwf = (
        ecmwf_store.select_time(t1, t2, n_levels=ecmwf_levels, combine='by_coords')
        .mean(dim='time')
        .swap_dims({"level":"pressure"})
    )

    ds_ecmwf = read_add_geopotential_altitude(ds_ecmwf)
    # print('ECMWF min pressure: ', min(ds_ecmwf.pressure.values), ', corresponding to geometric_height = ', ds_ecmwf.geometric_height.sel(pressure=min(ds_ecmwf.pressure.values)).values)
    #  print('ECMWF max pressure: ', max(ds_ecmwf.pressure.values), ', corresponding to geometric_height = ', ds_ecmwf.geometric_height.sel(pressure=max(ds_ecmwf.pressure.values)).values)

    return ds_ecmwf

def read_o3_apriori_ecmwf_mls_gromosOG_old(filename):
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

def read_o3_apriori_ecmwf_mls_gromosOG(folder_name, month_num):
    '''
    read the apriori o3 used in gromos retrieval
    
    Combined ECMWF and MLS climatology, described in:
    Studer et al. 2013, AMTD
    '''
    filename = folder_name + 'apriori_gromos_new_' + str(month_num) + '.aa'
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

    o3_apriori_ds = xr.Dataset({'o3': (['altitude','month'], 1e-6*ds.iloc[:,1:])},
                            coords = {
                                'altitude' : ('altitude', 1e3*ds.iloc[:,0]),
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
    month = pd.to_datetime(retrieval_param['time']).month

    if month in summer_months:
        fascod_clim = 'midlatitude-summer'
    else:
        fascod_clim = 'midlatitude-winter'
    
    fascod_atm = arts.Atmosphere.from_arts_xml(
        ARTS_DATA_PATH + "/planets/Earth/Fascod/{}/{}.".format(fascod_clim,fascod_clim)
    )
    print('Atmospheric state and apriori defined from Fascod climatology : ',fascod_clim)

    if retrieval_param['atm'] == 'fascod_cira86':
        # Reading CIRA86
        cira86 = read_cira86_monthly(retrieval_param['cira86_path'], month, retrieval_param['lat'])
        fascod_atm.set_t_field(cira86['Pressure'].values, cira86['temperature'].values)
        # z field
        fascod_atm.set_z_field(cira86['Pressure'].values, cira86["altitude"].values)
        print('Set altitude and Temperature fields to CIRA86 !')
        
    if retrieval_param['o3_apriori']=='gromos':
        o3_apriori = read_o3_apriori_ecmwf_mls_gromosOG(retrieval_param['apriori_ozone_climatology_GROMOS'], month)
        fascod_atm.set_vmr_field(
            "o3", o3_apriori["p"].values, o3_apriori['o3'].values
        )
        print('Ozone apriori from : OG GROMOS')
    elif retrieval_param['o3_apriori']=='somora':
        # For SOMORA, we need to interpolate the o3 apriori with altitude from FM (no pressure defined)
        o3_apriori = read_o3_apriori_OG_SOMORA(retrieval_param['apriori_ozone_climatology_SOMORA'], month)
        o3_apriori_h = interpolate(
            fascod_atm.z_field.data[:,0,0], 
            o3_apriori.altitude.data,
            o3_apriori.o3.data
            )
        # extracting pressure from the fascod atm
        pressure_fascod = fascod_atm.z_field.to_xarray().coords['Pressure']
        fascod_atm.set_vmr_field(
            "o3", pressure_fascod, o3_apriori_h
        )
        print('Ozone apriori from : OG SOMORA')
    else:
        print('Ozone apriori from : fascod')
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

def read_merra2_data(retrieval_param, merra_location, t1, t2):
    '''
    Function to read MERRA2 monthly data
    
    Merra2 data reminder because of its weird dimensions names:

    phony_dim_0 = time
    phony_dim_1 = altitude
    phony_dim_2 = latitude ?
    phony_dim_3 = longitude ?
    phony_dim_4 = ??

    '''
    merra_location = '/mnt/tub/atmosphere/MERRA2/'
    t1 = datetime.date(2017,1,6)
    #merra_global_filename = merra_location + t1.strftime('%Y')+ '_' + t1.strftime('%m') +'/' + 'MERRA2_400.inst3_3d_asm_Nv.'+t1.strftime('%Y%m%d')+'.nc4'
    merra_BRN_filename = merra_location + 'BRN/'+ 'MERRA2_BRN_'+t1.strftime('%Y_%m')+'_diagnostic.h5'

    # merra2_monthly = xr.open_dataset(
    #     merra_monthly_filename,
    #     mask_and_scale=True,
    #     decode_times=True,
    #     decode_coords=True,
    #     #use_cftime=True,
    #     )

    #latitude = merra2_monthly.lat.sel(phony_dim_0=1, phony_dim_1=1)
    #ozone = merra2_monthly.O3.sel(phony_dim_0=1, phony_dim_2=65, phony_dim_3=65)
    #ozone.plot()

    merra2_monthly_info = xr.open_dataset(
        merra_BRN_filename,
        group='info',
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        #use_cftime=True,
        )

    merra2_monthly_variables = xr.open_dataset(
        merra_BRN_filename,
        group='wind',
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        #use_cftime=True,
        )

def get_apriori_atmosphere_fascod_ecmwf_cira86(retrieval_param, ecmwf_store, cira86_path, t1, t2, extra_time_ecmwf, z_grid=None):
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

    ecmwf_time1 = t1 - pd.Timedelta(extra_time_ecmwf, "h")
    ecmwf_time2 = t2 + pd.Timedelta(extra_time_ecmwf, "h")

    if retrieval_param['verbose'] > 1:
        print('Searching ECMWF data between: '+str(ecmwf_time1)+ 'and '+str(ecmwf_time2))

    # Read EACMWF data (oper for now)
    ECMWF_store_path = ecmwf_store
    ecmwf_prefix = retrieval_param['ecmwf_prefix']
    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix, ecmwf_time1, ecmwf_time2)

    if sum(np.isnan(ds_ecmwf.pressure.values)) > 1:
        raise ValueError('no ECMWF data')

    # Reading CIRA86
    cira86 = read_cira86_monthly(cira86_path, month, lat)

    #plot_ecmwf_cira86_profile(ds_ecmwf, cira86)
        
    #filename = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'

    # Merging ecmwf and CIRA86
    ds_ptz = merge_ecmwf_cira86(ds_ecmwf, cira86, retrieval_param, method=retrieval_param['ptz_merge_method'], max_T_diff=retrieval_param['ptz_merge_max_Tdiff'])

    # extrapolate to p_grid ??
    if z_grid is not None:
        ds_ptz = extrapolate_down_ptz(ds_ptz, z_grid)

    if retrieval_param['verbose'] > 1:
        print('Min pressure of the new grid: ', min(ds_ptz.p.values))
        print('Max pressure of the new grid: ', max(ds_ptz.p.values))

        print('Min altitude of the new grid (interpolated from cira86 !!): ', min(ds_ptz.z.values))
        print('Max altitude of the new grid (interpolated from cira86 !!): ', max(ds_ptz.z.values))
    
        plot_apriori_ptz(ds_ptz)

    #  Temperature
    atm.set_t_field(ds_ptz['p'].values, ds_ptz['t'].values)
    
    if retrieval_param["retrieval_type"] == 8:
        if retrieval_param["test_type"] == 'Tprofile':
            print('Adding bias on Tprofile !')
            atm.set_t_field(ds_ptz['p'].values, ds_ptz['t'].values+retrieval_param['Tprofile_bias'])


    # z field
    atm.set_z_field(ds_ptz["p"].values, ds_ptz["z"].values)
    
    pressure_atm = atm.z_field.to_xarray().coords['Pressure']
    # DO NOT ADD O3 from ECMWF --> no value over 2 Pa...
    # Ozone
    #print('Pressure [Pa]:')
    #print(pressure_atm)

    if retrieval_param['o3_apriori'] == 'somora':
        retrieval_param['apriori_ozone_climatology_SOMORA'] = '/storage/tub/instruments/gromos/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'
        o3_apriori_SOMORA = read_o3_apriori_OG_SOMORA(retrieval_param['apriori_ozone_climatology_SOMORA'], month)
        o3_apriori_h = interpolate(
            atm.z_field.data[:,0,0], 
            o3_apriori_SOMORA.altitude.data,
            o3_apriori_SOMORA.o3.data
        )
        print('Ozone apriori from : old SOMORA retrievals')
        # extracting pressure from the fascod atm
        atm.set_vmr_field(
            "O3", pressure_atm, o3_apriori_h
        )       
    elif retrieval_param['o3_apriori'] == 'gromos':
        retrieval_param['apriori_ozone_climatology_GROMOS'] = '/storage/tub/instruments/gromos/InputsRetrievals/apriori_ECMWF_MLS/'
        o3_apriori_GROMOS = read_o3_apriori_ecmwf_mls_gromosOG(retrieval_param['apriori_ozone_climatology_GROMOS'], month)
        print('Ozone apriori from : old GROMOS retrievals')
        atm.set_vmr_field(
            "O3", o3_apriori_GROMOS["p"].values, o3_apriori_GROMOS['o3'].values
        )
    elif retrieval_param['o3_apriori'] == 'waccm':
        print('Ozone apriori from : WACCM climatology')
        ds_waccm = read_waccm(retrieval_param)
        atm.set_vmr_field(
            "O3", ds_waccm["p"].values, ds_waccm['o3'].values
        )
    elif retrieval_param['o3_apriori'] == 'waccm_monthly':
        print('Ozone apriori from : WACCM monthly (day/night) climatology')
        ds_waccm = read_waccm_monthly(retrieval_param)
        atm.set_vmr_field(
            "O3", ds_waccm["p"].values, ds_waccm['o3'].values
        )
    elif retrieval_param['o3_apriori'] == 'mls':
        o3_apriori = read_mls(retrieval_param['o3_apriori_file'])
        print('Ozone apriori from : mean MLS profile')
        atm.set_vmr_field(
            "O3", o3_apriori["p"].values, o3_apriori['o3'].values
        )
    elif retrieval_param['o3_apriori'] == 'retrieved_gromos':
        o3_apriori = read_retrieved(retrieval_param['o3_apriori_file'])
        ind =  np.where((o3_apriori.o3_p<10000) & (o3_apriori.o3_p>2.54))
        o3_apriori = o3_apriori.isel(o3_p=ind[0])
        print('Ozone apriori from : mean GROMOS profile, cut for MLS pressure grid')
        atm.set_vmr_field(
            "O3", o3_apriori["o3_p"].values, o3_apriori['o3_x'].values
        )
    elif retrieval_param['o3_apriori'] == 'retrieved_somora':
        o3_apriori = read_retrieved(retrieval_param['o3_apriori_file'])
        ind =  np.where((o3_apriori.o3_p<10000) & (o3_apriori.o3_p>2.54))
        o3_apriori = o3_apriori.isel(o3_p=ind[0])
        print('Ozone apriori from : mean SOMORA profile, cut for MLS pressure grid')
        atm.set_vmr_field(
            "O3", o3_apriori["o3_p"].values, o3_apriori['o3_x'].values
        )
    else: 
        print('Ozone apriori not recognized !')
        return 0
   # compare_o3_apriori_OG(o3_apriori_GROMOS.p.data, o3_apriori_GROMOS.o3, pressure_atm.data, o3_apriori_h)

    # Water vapor
    # merging Fascod with ECMWF
    if retrieval_param['h2o_apriori']=='fascod_ecmwf':
        fascod_atm = arts.Atmosphere.from_arts_xml(
            ARTS_DATA_PATH + "/planets/Earth/Fascod/{}/{}.".format(fascod_clim,fascod_clim)
        )
        h2o_apriori = merge_ecmwf_Fascod_atm(ds_ecmwf, fascod_atm)
        atm.set_vmr_field(
        "H2O", pressure_atm, h2o_apriori
        )
    elif retrieval_param['h2o_apriori']=='ecmwf':
        #plot_apriori_ptz(h2o_apriori)
        atm.set_vmr_field(
        "H2O", ds_ecmwf["pressure"].values, ds_ecmwf['specific_humidity'].values
        )
        if retrieval_param['verbose'] > 1:
            print('h2o taken from ecmwf')
    elif retrieval_param['h2o_apriori']=='fascod_extended':
        fascod_atm = arts.Atmosphere.from_arts_xml(
            ARTS_DATA_PATH + "/planets/Earth/Fascod/{}/{}.".format(fascod_clim,fascod_clim)
        )
        fascod_h2o = fascod_atm.vmr_field('h2o').to_xarray()
        fascod_h2o = fascod_h2o.sel(Latitude=0, Longitude=0)
        h2o_apriori = p_interpolate(
            pressure_atm,
            fascod_atm.z_field.to_xarray().coords['Pressure'].data, 
            fascod_h2o.data,
            )  
        atm.set_vmr_field(
        "H2O", pressure_atm, h2o_apriori
        )
        if retrieval_param['verbose'] > 1:
            print('h2o taken from fascod (extended)')
    else:
        raise ValueError('select valid apriori for h2o')
        
    print('Atmospheric state defined with: ECMWF oper v2, CIRA86')

    return atm

def read_mls(filename):
    mls_o3 = xr.open_dataset(
        filename,
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        )

    mls_o3['p'] = mls_o3['p']*100
    mls_o3['o3'] = mls_o3['o3']*1e-6
    return mls_o3

def read_waccm(retrieval_param, extra_day=0):
    filename = retrieval_param['waccm_file']
    datetime = retrieval_param['time']

    ds_waccm = xr.open_dataset(
        filename,
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        )

    # Introduce the solar zenith angle to decide for the apriori:
    lst, ha, sza, night = get_LST_from_GROMORA(datetime, retrieval_param['lat'], retrieval_param['lon'])
    
    #(sza,day,night) = solar_zenith_angle(datetime,retrieval_param)
    if night:
        tod = 'night'
    else:
        tod = 'day' 
    print('Solar elevation angle = ',90-sza, ', using ',tod,'time apriori profile !')
    # As a function of datetime, select appropriate day and time from climatology:
    if extra_day:
        ds_waccm = ds_waccm.sel(
            time=slice(pd.to_datetime(datetime).dayofyear-extra_day,pd.to_datetime(datetime).dayofyear+extra_day), 
            tod=tod
        ).mean(dim='time')
    else:
        ds_waccm = ds_waccm.sel(time=pd.to_datetime(datetime).dayofyear, tod=tod)
    return ds_waccm

def solar_zenith_angle_old(time,retrieval_param):
    sec_after_midnight = pd.to_datetime(time).hour*3600 + pd.to_datetime(time).minute*60 + pd.to_datetime(time).second 
    sec_to_noon = sec_after_midnight - 12*3600
    cos_solar_hour_angle = np.cos(np.deg2rad((sec_to_noon/3600)*15))

    # sun declination
    declination = -23.44*np.cos(np.deg2rad((pd.to_datetime(time).dayofyear+10)*360/365))
    cos_declination = np.cos(np.deg2rad(declination))
    sin_declination = np.sin(np.deg2rad(declination))

    # Sunrise/Sunset:
    cos_hour_angle_night= -np.tan(np.deg2rad(retrieval_param['lat']))*np.tan(np.deg2rad(declination))
    if cos_solar_hour_angle < cos_hour_angle_night:
        night = True
    else:
        night = False
    cos_sza = sin_declination*np.sin(np.deg2rad(retrieval_param['lat'])) + np.cos(np.deg2rad(retrieval_param['lat']))*cos_declination*cos_solar_hour_angle
    sza = np.rad2deg(np.arccos(cos_sza))
    if sza > 90:
        day = False
    else:
        day = True

    date = dt.datetime(
        pd.to_datetime(time).year,
        pd.to_datetime(time).month,
        pd.to_datetime(time).day,
        pd.to_datetime(time).hour,
        pd.to_datetime(time).minute,
        pd.to_datetime(time).second,
        tzinfo=timezone('UTC')
        )
    
    #date = date.tz_localize('Europe/Zurich')
    # Using pysolar package:
    sza_pysolar = solar.get_altitude(retrieval_param['lat'], retrieval_param['lon'], date)
    print('solar elevation angle = ', sza_pysolar)
    return sza, day, night

def read_waccm_monthly(retrieval_param):
    filename = retrieval_param['waccm_file']
    datetime = retrieval_param['time']

    ds_waccm = xr.open_dataset(
        filename,
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        )

    # Introduce the solar zenith angle to decide for the apriori:
    lst, ha, sza, night, tc = get_LST_from_GROMORA(datetime, retrieval_param['lat'], retrieval_param['lon'])

    if night:
        tod = 'night'
    else:
        tod = 'day' 


    month = pd.to_datetime(datetime).month
    month_str = pd.to_datetime(datetime).strftime('%b')

    # As a function of datetime, select appropriate month from climatology:
    ds_waccm_monthly = ds_waccm.sel(
        time=slice(month_start[month-1],month_stop[month-1]), 
        tod=tod
    ).mean(dim='time')

    if retrieval_param['verbose'] > 0:
        print(f'WACCM {month_str} profile, sza = {sza:1f}, using {tod} time apriori profile !')

    return ds_waccm_monthly

def read_waccm_yearly(filename, datetime):
    ds_waccm = xr.open_dataset(
        filename,
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        )

    tod = 'day' 
    
    # As a function of datetime, select appropriate month from climatology:
    ds_waccm_yearly = ds_waccm.sel(
        tod=tod
    ).mean(dim='time')

    return ds_waccm_yearly
def read_retrieved(filename):
    retrieved_o3 = xr.open_dataset(
        filename,
        mask_and_scale=True,
        decode_times=True,
        decode_coords=True,
        )

    retrieved_o3['o3_p'] = retrieved_o3['o3_p']*100
    retrieved_o3['o3_x'] = retrieved_o3['o3_x']*1e-6
    return retrieved_o3

def extrapolate_down_ptz(ds_ptz, z_grid):
    '''
    extrapolate to get low enough

    '''
    new_p = interpolate(z_grid, ds_ptz.z.values, ds_ptz.p.values)
    new_t = interpolate(z_grid, ds_ptz.z.values, ds_ptz.t.values)

    ds_extrapolated = xr.Dataset({'t': ('p', new_t),
                        'z': ('p', z_grid)},
                        coords = {
                            'p' : ('p', new_p),
                        }
    )
    return ds_extrapolated

def merge_ecmwf_cira86(ds_ecmwf, cira86, retrieval_param, method='simple_stack_corr', max_T_diff=10, max_pressure=50):
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

    if method=='simple':
        #upper_p_grid = cira86.Pressure.data[cira86.Pressure < np.min(ds_ecmwf.pressure.data)

        # selection based on altitude
        upper_p_grid = cira86.Pressure.data[cira86.altitude > np.max(ds_ecmwf.geometric_height.data)]
        #upper_t = cira86.temperature.data[cira86.Pressure < np.min(ds_ecmwf.pressure.data)]

        upper_cira86_ds = cira86.sel(Pressure = upper_p_grid)
        p_grid = np.hstack((ds_ecmwf.pressure.data, upper_p_grid))
        temperature = np.hstack((ds_ecmwf.temperature.data, upper_cira86_ds.temperature.data))
        z_grid = np.hstack((ds_ecmwf.geometric_height.data, upper_cira86_ds.altitude.data))
    elif method=='simple_stack_corr':
        # selection based on altitude
        upper_cira86_ds = cira86.where(cira86.Pressure < np.min(ds_ecmwf.pressure.data), drop=True)
        #upper_t = cira86.temperature.data[cira86.Pressure < np.min(ds_ecmwf.pressure.data)]

        #upper_cira86_ds = cira86.sel(Pressure = upper_p_grid)
        #ecmwf_lower = ds_ecmwf.where(ds_ecmwf.pressure > p_with_low_T_diff[-1])

        p_grid = np.hstack((ds_ecmwf.pressure.data, upper_cira86_ds.Pressure))
        temperature = np.hstack((ds_ecmwf.temperature.data, upper_cira86_ds.temperature.data))
        z_grid = np.hstack((ds_ecmwf.geometric_height.data, upper_cira86_ds.altitude.data))

    elif method=='max_diff':
        ecmwf_t_i = p_interpolate(
            cira86.Pressure.data, ds_ecmwf.pressure.data, ds_ecmwf.temperature.data, fill=np.nan
        )
        T_diff = cira86.temperature-ecmwf_t_i
    
        p_with_low_T_diff = T_diff.Pressure.where((T_diff<max_T_diff), drop=True)
        upper_p_grid = cira86.Pressure.where(cira86.Pressure<p_with_low_T_diff[-1], drop=True)
        upper_cira86_ds = cira86.sel(Pressure = upper_p_grid)
        ecmwf_lower = ds_ecmwf.where(ds_ecmwf.pressure > p_with_low_T_diff[-1])
        # # interpolate CIRA86 altitude on p_grid
        # cira86_alt_i = p_interpolate(
        #     p_grid, cira86.Pressure.data, cira86.altitude.data, fill = np.nan
        # )
        p_grid = np.hstack((ecmwf_lower.pressure.data, upper_p_grid))
        temperature = np.hstack((ecmwf_lower.temperature.data, upper_cira86_ds.temperature.data))
        z_grid = np.hstack((ecmwf_lower.geometric_height.data, upper_cira86_ds.altitude.data))
    elif method=='max_diff_surf':
        ecmwf_t_i = p_interpolate(
            cira86.Pressure.data, ds_ecmwf.pressure.data, ds_ecmwf.temperature.data, fill=np.nan
        )
        T_diff = cira86.temperature-ecmwf_t_i
    
        p_with_low_T_diff = T_diff.Pressure.where((T_diff<max_T_diff), drop=True)
        upper_p_grid = cira86.Pressure.where(cira86.Pressure<p_with_low_T_diff[-1], drop=True)
        upper_cira86_ds = cira86.sel(Pressure = upper_p_grid)
        ecmwf_lower = ds_ecmwf.where(ds_ecmwf.pressure > p_with_low_T_diff[-1])
        # # interpolate CIRA86 altitude on p_grid
        # cira86_alt_i = p_interpolate(
        #     p_grid, cira86.Pressure.data, cira86.altitude.data, fill = np.nan
        # )
        p_grid = np.hstack((retrieval_param['p_surface']*100, ecmwf_lower.pressure.data, upper_p_grid))
        temperature = np.hstack((retrieval_param['T_surface'], ecmwf_lower.temperature.data, upper_cira86_ds.temperature.data))
        z_grid = np.hstack((retrieval_param["station_altitude"], ecmwf_lower.geometric_height.data, upper_cira86_ds.altitude.data))

    ds_merged = xr.Dataset({'t': ('p', temperature),
                            'z': ('p', z_grid)},
                            coords = {
                                'p' : ('p', p_grid),
                            }
    )

    return ds_merged

def get_ph_levs(sp, level):
    '''Return the presure at a given level and the next'''
    #a_coef = values['pv'][0:values['nlevels'] + 1]
    #b_coef = values['pv'][values['nlevels'] + 1:]

    ph_lev = levels.hybrid_level_a[level-2] + (levels.hybrid_level_b[level-2] * sp)
    ph_levplusone = levels.hybrid_level_a[level-1] + (levels.hybrid_level_b[level-1] * sp)
    return ph_lev, ph_levplusone

def compute_z_level(lev, ds_ecmwf, z_h):
    '''Compute z at half & full level for the given level, based on t/q/sp'''
    # select the levelist and retrieve the vaules of t and q
    # t_level: values for t
    # q_level: values for q
    # codes_index_select(idx, 'level', lev)
    # codes_index_select(idx, 'shortName', 't')
    # gid = codes_new_from_index(idx)
    # if gid is None:
    #     raise MissingLevelError('T at level {} missing from input'.format(lev))
    # t_level = codes_get_values(gid)
    # codes_release(gid)
    # codes_index_select(idx, 'shortName', 'q')
    # gid = codes_new_from_index(idx)
    # if gid is None:
    #     raise MissingLevelError('Q at level {} missing from input'.format(lev))
    # q_level = codes_get_values(gid)
    # codes_release(gid)
    R_D = 287.06
    R_D = 287.0597
    R_G = 9.80665
    level_ind = np.where(ds_ecmwf.level.values == lev)[0][0]
    
    t_level = ds_ecmwf.temperature[level_ind].values
    q_level = ds_ecmwf.specific_humidity[level_ind].values
    
    # compute moist temperature
    t_level = t_level * (1. + 0.609133 * q_level)

    sp = np.exp(ds_ecmwf.logarithm_of_surface_pressure.values) 

    # compute the pressures (on half-levels)
    ph_lev, ph_levplusone = get_ph_levs(sp, lev)
 
    if lev == 1:
        dlog_p = np.log(ph_levplusone / 0.1)
        alpha = np.log(2)
    else:
        dlog_p = np.log(ph_levplusone / ph_lev)
        alpha = 1. - ((ph_lev / (ph_levplusone - ph_lev)) * dlog_p)

    t_level = t_level * R_D
 
    # z_f is the geopotential of this full level
    # integrate from previous (lower) half-level z_h to the
    # full level
    z_f = z_h + (t_level * alpha)
 
    # z_h is the geopotential of 'half-levels'
    # integrate z_h to next half level
    z_h = z_h + (t_level * dlog_p)
 
    return z_h, z_f
  
def read_add_geopotential_altitude(ds_ecmwf):
    '''Compute z at half & full level for the given level, based on t/q/sp'''
    # We want to integrate up into the atmosphere, starting at the
    # ground so we start at the lowest level (highest number) and
    # keep accumulating the height as we go.
    # See the IFS documentation, part III
    # For speed and file I/O, we perform the computations with
    # numpy vectors instead of fieldsets.
 
    z_h = ds_ecmwf.geopotential.values
    #codes_set(values['sample'], 'step', int(step))
    geopotential = -999*np.ones(len(ds_ecmwf.pressure))
    geometric_height = -999* np.ones(len(ds_ecmwf.pressure))
    i = -1
    for lev in sorted(ds_ecmwf.level.data, reverse=True):
        i=i+1
        #print(i)
        try:
            z_h, z_f = compute_z_level(lev, ds_ecmwf, z_h)
            geopotential[i]=z_f
            geometric_height[i]=z_f/9.80665
            # store the result (z_f) in a field and add to the output
            #codes_set(values['sample'], 'level', lev)
            #codes_set_values(values['sample'], z_f)
            #codes_write(values['sample'], fout)
        except MissingLevelError as e:
            print('%s [WARN] %s' % (sys.argv[0], e),
                  file=sys.stderr)


    ds = xr.Dataset(
            {
            'geometric_height': ('pressure', geometric_height.data),
            'geopotential_ecmwf': ('pressure', geopotential.data)
            },
            coords = {
            'pressure' : ('pressure', ds_ecmwf.pressure.data),
            }
    )

    merged = ds_ecmwf.merge(ds)
    return merged
    
def ecmwf_zp_calc( ds_ecmwf):
    '''
     [ z, p ] = ecmwf_zp_calc( A_h, B_h, T, p_surf, Phi_surf, q )
    
     Calculates geometric altitude and pressure from ECMWF profiles on model
     levels.
    
     Input parameters:
    
     A_h:      A_k+1/2 constant from GDS section of GRIB file on half levels
     B_h:      B_k+1/2 constant from GDS section of GRIB file on half levels
     T:        Temperature [K] on full levels
     p_surf:   Surface pressure [Pa] (optional)
     Phi_surf: Surface geopotential [m^2 s^-2] (optional)
     q:        Specific humidity [kg kg^-1] on full levels (optional)
    
     If missing, the parameters p_surf, Phi_surf, and q will be initialized
     to p_surf = 1.01325e5 Pa, Phi_surf = 0 m^2/s^2, q = 0 kg/kg.
    
     Output parameters:
    
     z: geometric altitude [m] on full levels
     p: pressure [Pa] on full levels
    
     Note: the function assumes that all vectors are sorted so that the first
           value is near the top of the atmosphere and the last value is near
           the surface.
    
     $Id: ecmwf_zp_calc.m,v 1.1 2005/03/04 23:06:53 feist Exp $
    
     Nomenclature for full and half-level variables:
    
     X is a variable on full levels (like all the prognostic variables)
     X_h is a variable on half levels
    
     In the ECMWF documentation in chapter 2.2 that would be written as
     X   <-> X_k
     X_h <-> X_k+1/2 (which is not possible in Matlab syntax)
    '''
    # # Set missing variables to default values

    # # if ~exist('p_surf', 'var')
    # #   p_surf = 100*1013.250; % Pa
    # # end

    # # if ~exist('Phi_surf', 'var')
    # #   Phi_surf = 0;
    # # end

    # # if ~exist('q', 'var')
    # #   q = zeros(size(T)); % Specific humidity [kg/kg]
    # # end

    # #
    # # Define physical constants according to the subroutine SUCST.F90
    # # from the source code of the ECMWF ITS model
    # #
    # RKBOL = 1.380658E-23;	# Boltzmann's constant k [J/K]
    # RNAVO = 6.0221367E+23;	# Avogadro's number NA []
    # R = RNAVO * RKBOL;	# ideal gas constant R [J/(mol*K)]
    # RMD = 28.9644;		# Dry air molecular weight [g/mol]
    # RMV = 18.0153;          # Water vapor molecular weight [g/mol]
    # RD = 1000 * R / RMD;	# Dry air constant Rd [J/(K*kg)]
    # RV = 1000 * R / RMV;    # Water vapor constand Rv [J/(K*kg)]
    # RG = 9.80665;		# Earth's gravitational acceleration g [m/s^2]

    # # # Get number of levels
    # NLEV = len(ds_ecmwf.pressure)

    
    # #Calculate half-level pressure p_h  from
    # #A_h, B_h and surface pressure p_surf
    # #Level 1 = top, NLEV = bottom
    # p_surf = np.exp(ds_ecmwf.logarithm_of_surface_pressure.data)
    # p_h = levels.hybrid_level_a + levels.hybrid_level_b * p_surf;              # Half-level pressure (eq. 2.11)
    # #p_h2 = ds_ecmwf.pressure.values
    # #
    # # Calculate delta_p according to Eq. 2.13
    # #
    # #delta_p = np.diff( p_h )

    # #
    # # Calculate virtual temperature according to standard textbook formula
    # #
    # q = np.flip(ds_ecmwf.specific_humidity.values)
    # T = np.flip(ds_ecmwf.temperature.values)
    # T_v = T*( 1 + ( RV / RD - 1 ) * q ); # Virtual temperature [K]

    # phi_h = ds_ecmwf.geopotential.data

    # # reversedLevel = np.flip(np.arange(0,137))
    # # for rlev in reversedLevel:
    # #     if rlev == max(reversedLevel):

    # #     else:
        
    # #     phi_h_new = phi_h + 
        

    # #
    # # Calculate ln_p = log( p_k+1/2 / p_k-1/2 ) for 1 <= k <= NLEV
    # #
    # ln_p = np.log( p_h[1:-1] / p_h[0:-2])
    # ln_p[0] = 0; # ln_p(1) is never used and might be infinite

    # #
    # # Geopotential height at half levels (eq. 2.21 in matrix form)
    # #
    # loopj = triu( np.ones( NLEV+1, NLEV ));  # Loop over j=k+1...NLEV as a matrix
    # Phi_h = Phi_surf + RD * loopj * ( T_v * ln_p );
    # Phi_h(1) = NaN; # Phi_h(1) was calculated with ln_p(1) and is never used

    # #
    # # Calculate alpha according to eq. 2.23
    # #
    # alpha = 1 - p_h[0:-1]/ delta_p * ln_p
    # alpha[0] = np.log(2); # Top level specially defined

    # #
    # # Calculate geopotential on full levels according to eq. 2.22
    # #
    # Phi = Phi_h(2:end) + RD * alpha .* T_v;

    # #
    # # Return geometric height and pressure on full levels
    # #
    # p = ( p_h(1:end-1) + p_h(2:end) ) / 2; # Full level pressure [Pa]
    # z = Phi / RG; # Geometric height [m]

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

    fig.suptitle('Apriori raw ptz profile')

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

def compare_o3_apriori_OG(pressure, o3_apriori, pressure_SOMORA, o3_apriori_SOMORA):
    fig, axs = plt.subplots(1,2, sharex=True)
    axs[0].plot(o3_apriori*1e6, pressure)
    axs[0].invert_yaxis()
    axs[0].set_yscale('log')
    axs[0].set_xlabel('$O_3$ [PPM]')
    axs[0].set_ylabel('$P$ [Pa]')
    axs[0].set_ylim(1e5,1e-3)
    axs[0].set_title('GROMOS')

    axs[1].plot(o3_apriori_SOMORA*1e6, pressure_SOMORA)
    axs[1].invert_yaxis()
    axs[1].set_yscale('log')
    axs[1].set_xlabel('$O_3$ [PPM]')
    axs[1].set_ylabel('$P$ [Pa]')
    axs[1].set_ylim(1e5,1e-3)
    axs[1].set_title('SOMORA')

    fig.suptitle('OG $O_3$ apriori')

    fig.show()

def read_ptz_ecmwf_cira86(time, lat, location, ecmwf_store, cira86_path, merge_method, max_T_diff, extra_time_ecmwf = 3.5, z_grid=None):
    '''
    Reading the a-priori atmosphere defined from ECMWF operationnal dataset and CIRA86 clim

    This function is a work around to add temperature profile information on the retrieved profiles.

    time: datetime
    '''
    summer_months = [4, 5, 6, 7, 8, 9]


    month = pd.to_datetime(time).month

    # Range of time to seach for ECMWF profiles.
    ecmwf_time1 = time - pd.Timedelta(extra_time_ecmwf, "h")
    ecmwf_time2 = time + pd.Timedelta(extra_time_ecmwf, "h")

    # Read EACMWF data (oper for now)
    ECMWF_store_path = ecmwf_store
    ecmwf_prefix = ecmwf_prefix = f'ecmwf_oper_v{2}_{location}_%Y%m%d.nc'
    ds_ecmwf = extract_ecmwf_ds(ECMWF_store_path, ecmwf_prefix, ecmwf_time1, ecmwf_time2)

    if sum(np.isnan(ds_ecmwf.pressure.values)) > 1:
        raise ValueError('no ECMWF data')

    # Reading CIRA86
    cira86 = read_cira86_monthly(cira86_path, month, lat)

    #plot_ecmwf_cira86_profile(ds_ecmwf, cira86)
        
    #filename = '/home/eric/Documents/PhD/GROSOM/InputsRetrievals/AP_ML_CLIMATO_SOMORA.csv'

    # Merging ecmwf and CIRA86
    ds_ptz = merge_ecmwf_cira86(ds_ecmwf, cira86, method=merge_method, max_T_diff=max_T_diff)
    

    
    # extrapolate to p_grid ??
    if z_grid is not None:
        ds_ptz = extrapolate_down_ptz(ds_ptz, z_grid)

    return ds_ptz
    
if __name__ == "__main__":
    date = dt.datetime(2020,1,1)
    ptz = read_ptz_ecmwf_cira86(
        time=pd.to_datetime(date), 
        lat=46, 
        location='BERN', 
        ecmwf_store='/storage/tub/atmosphere/ecmwf/locations/'+str(date.year), 
        cira86_path=os.path.join(ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly'), 
        merge_method='max_diff',
        max_T_diff=5,
        extra_time_ecmwf = 3.5, 
        z_grid=None)