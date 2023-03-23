#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 

@author: eric

Script that contains the function to convert GROMORA level 2 to NDACC GEOMS compliant file


Including : 
    It includes some workaround function copied from gromora_atmosphere in order to complete the missing data in the current level 2 (like the Tprofile...)
"""
import sys, os, time
from os.path import dirname, abspath, join
# Adding the paths to pyretrievals and retrieval folder
sys.path.append(join(dirname(sys.path[0]),'pyretrievals'))
sys.path.append(join(dirname(sys.path[0]),'retrieval'))
from abc import ABC
import datetime as dt
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from scipy.interpolate import interp2d
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
from pytz import timezone
from gromora_utils import *
from gromora_time import *

import math

from pyhdf.SD import SD,SDC
import fnmatch,shutil,glob,scipy
from ndacc_parameters import *

# # For ARTS, we need to specify some paths
load_dotenv('/opt/arts/.env.stockhorn-arts24')
# load_dotenv('/opt/anaconda/.env.birg-arts24')
# load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')

#ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
#ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
#ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

from retrievals import arts
from retrievals.data.ecmwf import ECMWFLocationFileStore
from retrievals.data.ecmwf import levels
from retrievals.data import interpolate
from retrievals.data import p_interpolate
from retrievals.arts.atmosphere import p2z_simple, z2p_simple
#from typhon.arts.xml import load
from pyarts.xml import load
from GROMORA_library import cmap

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

def merge_ecmwf_cira86(ds_ecmwf, cira86, method='simple_stack_corr', max_T_diff=10, max_pressure=50):
    '''
    Merging profile from ECMWF oper and CIRA86 monthly climatology.

    Start with a very simple scheme, we just take CIRA86 values when we reach the top
    of the ECMWF ones.
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
        ecmwf_lower = ds_ecmwf.where(ds_ecmwf.pressure >= p_with_low_T_diff[-1], drop=True)
        # # interpolate CIRA86 altitude on p_grid
        # cira86_alt_i = p_interpolate(
        #     p_grid, cira86.Pressure.data, cira86.altitude.data, fill = np.nan
        # )
        p_grid = np.hstack((ecmwf_lower.pressure.data, upper_p_grid))
        temperature = np.hstack((ecmwf_lower.temperature.data, upper_cira86_ds.temperature.data))
        z_grid = np.hstack((ecmwf_lower.geometric_height.data, upper_cira86_ds.altitude.data))


    ds_merged = xr.Dataset({'t': ('p', temperature),
                            'z': ('p', z_grid)},
                            coords = {
                                'p' : ('p', p_grid),
                            }
    )

    #ds_merged = ds_merged.interpolate_na(dim='t', method='linear', max_gap=3)

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
        except e:
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

def read_cira86_monthly(cira86_datapath, month = None, latitude = None):
    '''
    Extract profile from the CIRA86 monthly climatology.

    For now, works with the ones saved in arts-data saved as .xml files:
    2 files per month in the form:
    cira86_month*.t.xml -> temperature data, zonally averaged
    cira86_month*.z.xml -> geometric altitude
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


def read_ptz_ecmwf_cira86(time, lat, location, ecmwf_store, cira86_path, merge_method, max_T_diff, extra_time_ecmwf = 3.5, z_grid=None):
    '''
    Reading the a-priori atmosphere defined from ECMWF operationnal dataset and CIRA86 clim

    This function is a work around to add temperature profile information on the retrieved profiles for the NDACC compliance.

    time: datetime
    '''
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

    # Merging ecmwf and CIRA86
    ds_ptz = merge_ecmwf_cira86(ds_ecmwf, cira86, method=merge_method, max_T_diff=max_T_diff)
    
    # extrapolate to p_grid ??
    if z_grid is not None:
        ds_ptz = extrapolate_down_ptz(ds_ptz, z_grid)

    return ds_ptz

def gromora_level2_GEOMS(instrument_name= "GROMOS", date= dt.date(2021, 6 , 27), spectros = ['AC240'] , ex = '_v2', new_z=1e3*np.arange(10, 90, 2), avk_corr = True, plot_tprofile=False, save_nc=True, RD=True, outfolder = '/home/es19m597/Documents/GROMORA/NDACC/'):
    """
    Function to convert GROMORA level 2 to GEOMS compliant xarray dataset
    As a option, it can save the dataset as netCDF4 file.


    Args:
        instrument_name (str, optional): Name of the instrument. Defaults to "GROMOS".
        date (datetime, optional): the date to treat. Defaults to dt.date(2021, 6 , 27).
        plot_tprofile (bool, optional): a boolean to plot the comparisons between the og and the NDACC file. Defaults to False.
        save_nc (bool, optional): Boolean to save the dataset in netCDF format.
        outfolder (str): folder name to save the nc file

    Returns:
        xarray dataset: the new generated dataset, GEMOS compliant
    """
    integration_strategy = 'classic'
    
    int_time = 1
    new_L2 = True 

    measurement_response_limit = 0.75

    spectro = spectros[0]

    #############################################################################
    #############################################################################
    if instrument_name == "GROMOS":
        basename_lvl1 = os.path.join(
            '/storage/tub/instruments/gromos/level1/GROMORA/v2/', str(date.year))
        #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
        if new_L2:
            basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v2/"+str(date.year)
        else:
            basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v1/"+str(date.year)

        if spectro == 'AC240':
            import gromos_classes as gromos
    
            instrument = gromos.GROMOS_LvL2(
                date=date,
                basename_lvl1=basename_lvl1,
                basename_lvl2=basename_lvl2,
                integration_strategy=integration_strategy,
                integration_time=int_time
            )
        elif spectro == 'FB':
            import gromos_FB_classes as gromos
            instrument = gromos.GROMOS_FB_LvL2(
                date=date,
                basename_lvl1=basename_lvl1,
                basename_lvl2=basename_lvl2,
                integration_strategy=integration_strategy,
                integration_time=int_time
            )

        loc='Bern'
        
    elif instrument_name == "SOMORA":
        import somora_classes as somora
        basename_lvl1 = os.path.join(
            '/storage/tub/instruments/somora/level1/v2/', str(date.year))
        if new_L2:
        #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
            basename_lvl2 = "/storage/tub/instruments/somora/level2/v2/"+str(date.year)
        else:
            basename_lvl2 = "/storage/tub/instruments/somora/level2/v1/"+str(date.year)
        instrument = somora.SOMORA_LvL2(
            date=date,
            basename_lvl1=basename_lvl1,
            basename_lvl2=basename_lvl2,
            integration_strategy=integration_strategy,
            integration_time=int_time
        )
        loc='PAYERNE'

    #############################################################################
    # Read the original GROMORA level 2:
    level2_dataset = instrument.read_level2(
        spectrometers=spectros,
        extra_base=ex
    )

    # Define retrieval quality:
    good_data = xr.where((np.abs(level2_dataset[spectro].oem_diagnostics[:, 2] - 1) < instrument.cost_threshold(date.year) ) & (np.abs(level2_dataset[spectro].poly_fit_x[:,0].data)<instrument.polyfit_threshold), True, False)
    level2_dataset[spectro]['retrieval_quality'] = ('time', good_data.data.astype(int))

    if spectro == 'AC240':
        dataset = level2_dataset[spectro].drop_dims(['f', 'h2o_continuum_p', 'h2o_continuum_p_avk', 'poly_order', 'f_shift_grid', 'sine_grid'])
    elif spectro == 'FB':
        dataset = level2_dataset[spectro].drop_dims(['f', 'h2o_continuum_p', 'h2o_continuum_p_avk', 'poly_order', 'f_shift_grid'])

    #dataset['o3_p'] = dataset.o3_z.mean(dim='time').data
    #new_z = np.int32(dataset.o3_z.mean(dim='time').data)
    
    #.replace(tzinfo=gromora_tz)
    new_z = np.int32(new_z)

    lat = dataset.lat.mean(dim='time').data
    lon = dataset.lon.mean(dim='time').data  
    altitude_instrument = instrument.altitude 

    if not 'tropospheric_opacity' in list(dataset.data_vars.keys()) or not 'first_sky_time' in list(dataset.data_vars.keys()):
        # For the old level 2 format, we need to read the tropopheric opacity from the level 1b
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
            no_flag=False, meta_data=True, extra_base='')
        level1b = integrated_dataset[spectro]
        # Interpolating opactiy
        opacity = level1b.tropospheric_opacity.reindex_like(dataset.time, method='nearest', tolerance='1H').data
        first_time = level1b.first_sky_time.reindex_like(dataset.time, method='nearest', tolerance='1H').dt.round('s')
        last_time = level1b.last_sky_time.reindex_like(dataset.time, method='nearest', tolerance='1H').dt.round('s')
        number_of_spectra = level1b.number_of_sky_spectra.reindex_like(dataset.time, method='nearest', tolerance='1H').data
    else:
        opacity = dataset.tropospheric_opacity.data
        first_time = dataset.first_sky_time.dt.round('s')
        last_time = dataset.last_sky_time.dt.round('s')
        number_of_spectra = dataset.number_of_spectra.data

    #############################################################################
    # Saving the variables for each timestamps (easier for flagging single timestamp)

    o3 = list()
    o3_p = list()
    o3_random = list()
    o3_sytematic = list()
    o3_tot = list()
    o3_resolution = list()
    o3_apriori = list()
    o3_apriori_contribution = list()
    o3_avk = list()
    o3_number_density = list()
    o3_total_column =  list()
    temperature = list()
    mean_sza = list()

    for t in dataset.time.data:
        ds = dataset.sel(time=t)
        #if not 'solar_zenith_angle' in list(dataset.data_vars.keys()):
        #lst, ha, sza, night, tc = get_LST_from_GROMORA(datetime64_2_datetime(t).replace(tzinfo=gromora_tz),  lat, lon, check_format=False)
        lst, ha, sza, night, tc = utc2lst_NOAA(datetime64_2_datetime(t).replace(tzinfo=gromora_tz),  lat, lon)
        mean_sza.append(sza)
        #else:
        #    mean_sza.append(ds['solar_zenith_angle'].data)
        #print('DOF before interpolation = ', np.trace(ds.o3_avkm.data))
        # Interpolate all quantities on new altitude grid for NDACC
        ds.coords['alt'] = ('o3_p', ds.o3_z.data)
        
        ds = ds.swap_dims({'o3_p':'alt'})
        ds.coords['o3_p_avk'] = ('alt', ds.o3_z.data)
        
        #ds = ds.coords['o3_p_avk'].rename('alt')
        #ds=ds.drop_dims(['o3_p_avk'])
        
        #ds = ds.swap_dims({'o3_p_avk':'alt'})

        #f = np.interp(new_z, ds.o3_z.data[1:], np.diff(ds.o3_z.data))/np.median(np.diff(new_z))

        interp_ds = ds.interp(alt=new_z)

        # Specific 2D interpolation for the AVKs
        f_avkm_interp = interp2d(ds.o3_avkm.alt.data, ds.o3_avkm.o3_p_avk.data, ds.o3_avkm.data)
        avkm_interp = f_avkm_interp(new_z,new_z)

        # Correction factor for the AVKs interpolation
        # Initial sum of the AVKs and interpolation on new grid:
        slri = interp_ds.o3_mr.data #ds.o3_avkm.sum(axis=1).interp(alt=new_z) # Same as MR interpolated ! .where(ds.alt>new_z[0]).where(ds.alt<new_z[-1])

        shr = np.sum(avkm_interp, axis=1)
        factor =slri/shr

        avkm_interp_corr = np.ones_like(avkm_interp)*np.nan
        for h in range(len(avkm_interp)):
            avkm_interp_corr[h,:] = avkm_interp[h,:]*factor[h].data
            #avkm_interp_corr[h,:] = avkm_interp[h,:]*interp_ds.o3_mr.data[h]/np.sum(avkm_interp[h,:])

        # print('DOF after interpolation = ', np.trace(avkm_interp))
        # print('DOF difference = ', np.nansum(np.diag(ds.o3_avkm.where(ds.alt>new_z[0]).where(ds.alt<new_z[-1]).data))-np.trace(avkm_interp))
        # print('DOF difference (with corr) = ',np.nansum(np.diag(ds.o3_avkm.where(ds.alt>new_z[0]).where(ds.alt<new_z[-1]).data))-np.nansum(np.diag(avkm_interp_corr)))

        if not 'temperature_profile' in list(dataset.data_vars.keys()):
            ptz = read_ptz_ecmwf_cira86(
                time=pd.to_datetime(ds.time.data), 
                lat=instrument.latitude, 
                location=loc, 
                ecmwf_store='/storage/tub/atmosphere/ecmwf/locations/'+loc, 
                cira86_path=os.path.join(ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly'), 
                merge_method='max_diff',
                max_T_diff=5,
                extra_time_ecmwf = 3.5, 
                z_grid=new_z)
            #ds['temperature'] = (['time', 'o3_p'], ptz.t.data)

            # Interpolation of the dataset on the new altitude coordinates:
            t_interp = ptz.t.values #p_interpolate(ds.o3_p.data, ptz.p.values, ptz.t.values)
        else:
            t_interp = interp_ds['temperature_profile'].data

        # Clean the dataset for good quality data
        if ds.retrieval_quality.data == 1:
            ozone = interp_ds.o3_x.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data
            N_o3 = o3_vmr2number_density(ozone, interp_ds.o3_p,t_interp)
            o3_number_density.append(N_o3)
            tco=0
            for i in range(len(N_o3)-1):
                tco = tco + 0.5*(np.nansum([N_o3[i],N_o3[i+1]]))*(new_z[i+1]-new_z[i])
            o3_total_column.append(tco/DU)

            o3.append(1e6*ozone)
            o3_p.append(1e-2*interp_ds.o3_p.data)
            o3_random.append(1e6*interp_ds.o3_eo.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data)
            o3_sytematic.append(1e6*interp_ds.o3_es.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data)
            o3_tot.append(1e6*(interp_ds.o3_es+ interp_ds.o3_eo).where(interp_ds.o3_mr>measurement_response_limit, np.nan).data)

            o3_resolution.append((interp_ds.o3_fwhm.where(interp_ds.o3_mr>measurement_response_limit, -900000).data))
            o3_apriori_contribution.append(100*(1-interp_ds.o3_mr.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data))
        else:
            o3_number_density.append(np.ones_like(interp_ds.alt.data)*np.nan)
            o3_total_column.append(np.nan)
            o3.append(np.ones_like(interp_ds.alt.data)*np.nan)
            o3_p.append(np.ones_like(interp_ds.alt.data)*np.nan)
            o3_random.append(np.ones_like(interp_ds.alt.data)*np.nan)
            o3_sytematic.append(np.ones_like(interp_ds.alt.data)*np.nan)
            o3_tot.append(np.ones_like(interp_ds.alt.data)*np.nan)

            o3_resolution.append((interp_ds.o3_fwhm.where(interp_ds.o3_mr>measurement_response_limit, -900000).data))
            o3_apriori_contribution.append(100*(1-interp_ds.o3_mr.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data))
        if avk_corr:
            o3_avk.append(avkm_interp_corr)
        else:
            o3_avk.append(avkm_interp)
        o3_apriori.append(1e6*interp_ds.o3_xa)
        temperature.append(t_interp)
        #o3_mmr = o3_vmr2mmr(ozone)
    
    #interp_ds = dataset.interp(o3_p=z2p_simple(new_z))
    julian_dates = mjd2k_date(pd.to_datetime(dataset.time.data))

    start_date_iso = pd.to_datetime(
        first_time[0].data
        ).strftime('%Y%m%dT%H%M%SZ') #.round('1s')

    stop_date_iso = pd.to_datetime(
        last_time[-1].data
        ).strftime('%Y%m%dT%H%M%SZ') #.round('1s').
       
    first_time_MJD2K = np.float64(mjd2k_date(pd.to_datetime(first_time.data)))
    last_time_MJD2K = np.float64(mjd2k_date(pd.to_datetime(last_time.data)))

    # Check that the datetime are between first and last time:
    assert (julian_dates[0]>first_time_MJD2K[0]) and (julian_dates[-1]<last_time_MJD2K[-1])

    ####################################################
    # Add the GROMORA level 2 flags:
    flag=0
    for flagD in instrument.day2flag_level2:
        if len(dataset.time.sel(time=flagD))>0:     
            #print(flagD)      
            flag+=1
    if RD:
        data_quality = 'RD;DataNotQualityControlled'
    else:
        data_quality = 'Day is OK'
        if flag>0:
            data_quality = 'Day is not OK'
    # if opacity > 2:
    #     data_quality = data_quality + ', WARNING: high opacity'

    ####################################################
    # Global attributes definition
    file_version='001'
    #filename = 'groundbased_mwr.o3_'+instrument.affiliation+'_'+instrument.location.lower()+'_'+start_date_iso.lower()+'_'+stop_date_iso.lower()+'_'+file_version
    filename = 'groundbased_'+instrument.source.lower()+'_'+instrument.location.lower()+'_'+start_date_iso.lower()+'_'+stop_date_iso.lower()+'_'+file_version
    
    file_generation_date = pd.to_datetime(
        dt.datetime.now()
        ).round('1s').strftime('%Y%m%dT%H%M%SZ')

    # Global attributes dictionary:
    global_attrs = instrument.global_attributes_ndacc
    global_attrs['DATA_DISCIPLINE'] = 'ATMOSPHERIC.PHYSICS;REMOTE.SENSING;GROUNDBASED'
    global_attrs['DATA_GROUP'] = 'EXPERIMENTAL;PROFILE.STATIONARY'
    global_attrs['DATA_LOCATION'] = instrument.location
    global_attrs['DATA_SOURCE']=instrument.source
    global_attrs['DATA_START_DATE'] = start_date_iso
    global_attrs['DATA_STOP_DATE'] = stop_date_iso

    global_attrs['DATA_FILE_VERSION'] = file_version
    global_attrs['FILE_NAME']=filename+'.hdf'
    global_attrs['FILE_GENERATION_DATE'] = file_generation_date
    global_attrs['DATA_MODIFICATIONS'] = 'New harmonized retrievals from Swiss MWRs for FFTS time period described in: 10.5194/amt-15-6395-2022'
    global_attrs['DATA_CAVEATS'] = 'Ozone profiles are estimated with the Optimal Estimation Method as implemented in the Atmospheric Radiative Transfer Simulator (ARTS)'
    global_attrs['DATA_QUALITY'] = data_quality
    global_attrs['DATA_TEMPLATE'] = 'GEOMS-TE-MWR-003'
    global_attrs['DATA_PROCESSOR'] = 'Harmonized time series processed with '+dataset.arts_version
    global_attrs['FILE_ACCESS']='EVDC;NDACC'
    global_attrs['FILE_DOI']=' '
    global_attrs['FILE_ASSOCIATION']='NDACC'
    global_attrs['FILE_META_VERSION']='04R069;CUSTOM'

    ####################################################
    # Variables definition:
    # Written in a Python dictionary first for easy creation of the dataset later
    var_dict= dict(
        LATITUDE = (np.float32(interp_ds.lat.data)),
        LONGITUDE = (np.float32(interp_ds.lon.data)),
        ALTITUDE_INSTRUMENT = (np.int32(altitude_instrument)),
        ANGLE_VIEW_AZIMUTH = ('DATETIME', np.float32(dataset.obs_aa.data)),
        ANGLE_VIEW_ZENITH = ('DATETIME', np.float32(dataset.obs_za.data)),
        ANGLE_SOLAR_ZENITH_MEAN = ('DATETIME',np.float32(mean_sza)),
        #NOISE= ('DATETIME', dataset.median_noise.data),
        OPACITY_ATMOSPHERIC_EMISSION = ('DATETIME',np.float32(opacity)),
        DATETIME_START = ('DATETIME', first_time_MJD2K),
        DATETIME_STOP = ('DATETIME', last_time_MJD2K),
        INTEGRATION_TIME = ('DATETIME',np.float64(instrument.cycle_duration*number_of_spectra)),
        PRESSURE_INDEPENDENT = (['DATETIME', 'ALTITUDE'], np.float32(o3_p)),
        TEMPERATURE_INDEPENDENT = (['DATETIME', 'ALTITUDE'], np.float32(temperature)),
        O3_MIXING_RATIO_VOLUME_EMISSION = (['DATETIME', 'ALTITUDE'], np.float32(o3) ),
        O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_RANDOM_STANDARD = (['DATETIME', 'ALTITUDE'], np.float32(o3_random)),
        O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_SYSTEMATIC_STANDARD= (['DATETIME', 'ALTITUDE'], np.float32(o3_sytematic) ),
        O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_COMBINED_STANDARD = (['DATETIME', 'ALTITUDE'], np.float32(o3_tot) ),
        O3_MIXING_RATIO_VOLUME_EMISSION_RESOLUTION_ALTITUDE = (['DATETIME', 'ALTITUDE'], np.int32(o3_resolution)),
        O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI = (['DATETIME', 'ALTITUDE'], np.float32(o3_apriori) ),
        O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI_CONTRIBUTION = (['DATETIME', 'ALTITUDE'], np.float32(o3_apriori_contribution) ),
        O3_MIXING_RATIO_VOLUME_EMISSION_AVK = (['DATETIME', 'ALTITUDE', 'ALTITUDE'], np.float32(o3_avk) ),
        O3_COLUMN_PARTIAL_EMISSION = (['DATETIME'], np.float32(o3_total_column) ),
        O3_NUMBER_DENSITY_EMISSION = (['DATETIME', 'ALTITUDE'],  np.float32(o3_number_density) ),
        )

    ####################################################
    # Creation of the xarray dataset
    new_ds = xr.Dataset(
        coords=dict(
            DATETIME= julian_dates,
            ALTITUDE=('ALTITUDE', new_z),
        ),
        data_vars=var_dict,
        attrs=global_attrs
    )
    variables_all = ''
    varname = list(new_ds.keys())

    ####################################################
    # Variables attributes
    for s in list(new_ds.coords)+varname:
        n = list_indices[s]
        if s=='DATETIME':
            valid_min = np.float64(-2200)
            valid_max = np.float64(math.ceil(np.max(julian_dates)))
            fillV =  np.float64(-900000)
        elif s=='DATETIME_START':
            valid_min = np.float64(math.floor(np.min(first_time_MJD2K)))
            valid_max = np.float64(math.ceil(np.max(first_time_MJD2K)))
            fillV =  np.float64(-900000)
        elif s=='DATETIME_STOP':
            valid_min = np.float64(math.floor(np.min(last_time_MJD2K)))
            valid_max = np.float64(math.ceil(np.max(last_time_MJD2K)))
            fillV =  np.float64(-900000)
        else:
            if SDS_DataType_List[n] == 'REAL':
                valid_min = np.float32(VAR_VALID_MIN_List[n])
                valid_max = np.float32(VAR_VALID_MAX_List[n])
                fillV =  np.float32(-900000)
            elif SDS_DataType_List[n] == 'DOUBLE':
                valid_min = np.float64(VAR_VALID_MIN_List[n])
                valid_max = np.float64(VAR_VALID_MAX_List[n])
                fillV =  np.float64(-900000)
            elif SDS_DataType_List[n] == 'INTEGER':
                valid_min = np.int32(VAR_VALID_MIN_List[n])
                valid_max = np.int32(VAR_VALID_MAX_List[n])
                fillV = np.int32(-900000)

        if s =='O3_NUMBER_DENSITY_EMISSION':
            fillV = np.float32(-9e19)         

        if s == 'O3_MIXING_RATIO_VOLUME_EMISSION_AVK':
            avk = new_ds[s].data
            avk[np.isnan(new_ds[s].data)] = fillV
            avk[avk<valid_min] = fillV
            avk[avk>valid_max] = fillV
            new_ds[s] = (['DATETIME', 'ALTITUDE', 'ALTITUDE'], avk)
        # elif s == 'O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI':
        #     new_ds[s] = new_ds[s].where(~np.isinf(new_ds[s].data), fillV)
        #     new_ds[s] = new_ds[s].where(new_ds[s].data > valid_min, fillV)
        #     new_ds[s] = new_ds[s].where(new_ds[s].data < valid_max, fillV)
        #     new_ds[s] = new_ds[s].fillna(fillV) 
        else:
            new_ds[s] = new_ds[s].where(~np.isinf(new_ds[s].data), fillV)
            new_ds[s] = new_ds[s].where(new_ds[s].data > valid_min, fillV)
            new_ds[s] = new_ds[s].where(new_ds[s].data < valid_max, fillV)
            new_ds[s] = new_ds[s].fillna(fillV) 
        
        new_ds[s].attrs['VAR_NAME'] = SDS_Name_List[n]
        new_ds[s].attrs['VAR_DESCRIPTION'] = VAR_DESCRIPTION_List[n]
        new_ds[s].attrs['VAR_NOTES'] = VAR_NOTES_List[n]
        #new_ds[s].attrs['VAR_SIZE'] = VAR
        new_ds[s].attrs['VAR_DEPEND'] = VAR_DEPEND_List[n]
        new_ds[s].attrs['VAR_UNITS'] = VAR_UNITS_List[n]
        new_ds[s].attrs['VAR_SI_CONVERSION'] = VAR_SI_CONVERSION_List[n]
        new_ds[s].attrs['VAR_DATA_TYPE'] = SDS_DataType_List[n]
        new_ds[s].attrs['VAR_VALID_MIN'] = valid_min
        new_ds[s].attrs['VAR_VALID_MAX'] = valid_max
        new_ds[s].attrs['VAR_FILL_VALUE'] = fillV
        #new_ds[s].attrs['_FillValue'] = fillV
            
        if VAR_DEPEND_List[n] == 'CONSTANT':
            new_ds[s].attrs['VAR_SIZE'] = '1'
        elif VAR_DEPEND_List[n] == 'ALTITUDE':
            new_ds[s].attrs['VAR_SIZE'] = str(len(new_ds.ALTITUDE))
        elif VAR_DEPEND_List[n] == 'DATETIME':
            new_ds[s].attrs['VAR_SIZE'] = str(len(new_ds.DATETIME))
        elif VAR_DEPEND_List[n] == 'DATETIME;ALTITUDE':
            new_ds[s].attrs['VAR_SIZE'] = str(len(new_ds.DATETIME))+';'+str(len(new_ds.ALTITUDE))
        elif VAR_DEPEND_List[n] == 'DATETIME;ALTITUDE;ALTITUDE':
            new_ds[s].attrs['VAR_SIZE'] = str(len(new_ds.DATETIME))+';'+str(len(new_ds.ALTITUDE))+';'+str(len(new_ds.ALTITUDE))

        # rename variables according to NDACC conventions:
        new_ds = new_ds.rename_vars({s:SDS_Name_List[n]})
       
        variables_all = variables_all+SDS_Name_List[n]+';'

    variables_all = variables_all[:-1]
    new_ds.attrs['DATA_VARIABLES'] = variables_all

######################################################################################################################
# If True, save the dataset as netCDF files which should be close to GEOMS compliant. 
    if save_nc:
        new_ds.to_netcdf(outfolder+filename+'.nc')

######################################################################################################################
# Some plotting functions for checks    
    if plot_tprofile is not None:
        fig, axs = plt.subplots(nrows=1, ncols=4, sharey=True, figsize=(24,16))
        i = plot_tprofile[0]

        ds = level2_dataset[spectro].isel(time=i)
        new = new_ds.isel(DATETIME=i)
        ptz = read_ptz_ecmwf_cira86(
            time=pd.to_datetime(ds.time.data), 
            lat=instrument.latitude, 
            location=loc, 
            ecmwf_store='/storage/tub/atmosphere/ecmwf/locations/'+loc, 
            cira86_path=os.path.join(ARTS_DATA_PATH, 'planets/Earth/CIRA86/monthly'), 
            merge_method='max_diff',
            max_T_diff=5,
            extra_time_ecmwf = 3.5, 
            z_grid=None)

        # #p_grid = np.logspace(-3, 5, 47)
        # z_grid = p2z_simple(ds.o3_p.data)
        # p_grid = z2p_simple(np.arange(1e3, 95e3, 2e3))
        p_grid = ds.o3_p.data

        z_interp = p_interpolate(p_grid, ptz.p.values, ptz.z.values )
        t_interp = p_interpolate(p_grid, ptz.p.values, ptz.t.values)

        if not 'temperature_profile' in list(dataset.data_vars.keys()):
            axs[0].plot(ptz.t.values, ptz.z.values/1e3, 'x')
            axs[0].plot(t_interp, z_interp/1e3, 'o')
        else:
            axs[0].plot(ds.temperature_profile.data, ds.o3_z.values/1e3, 'x')
        axs[0].plot(new.TEMPERATURE_INDEPENDENT, new.ALTITUDE/1e3)
        axs[0].set_xlim(200,300)
        axs[0].set_ylim(0,95)

        o3_og = ds.o3_x.where(ds.o3_mr>measurement_response_limit, np.nan)
        o3 = new['O3.MIXING.RATIO.VOLUME_EMISSION']
        o3_z = ds.o3_z 
        new_alt = new.ALTITUDE
        o3_p = ds.o3_p
        y_lab = 'Altitude [km]'

        axs[3].plot(1e6*o3_og, o3_z/1e3, 'o-',label='OG')
        axs[3].plot(o3, new_alt/1e3,'x-', label='NDACC')
        axs[3].legend()
        axs[3].set_xlim(-0.1,10)
        #axs[1].fill_betweenx(y_axis, (o3-error)*to_ppm,(o3+error)*to_ppm, color=col, alpha=0.5)
        #axs[1].plot(o3*to_ppm, y_axis,'-', linewidth=1.5, label='retrieved',color=col)
        #axs[0].plot(o3_apriori*to_ppm, y_axis, '--', linewidth=1.5, label='apriori',color='k')
        #axs[0].set_title('O$_3$ VMR')
        #axs[0].set_xlim(-0.5,9)
        axs[1].plot(0.25*ds.o3_mr, o3_z/1e3, 'k-') 
        counter=0
        color_count = 0
        for j, avk in enumerate(ds.o3_avkm):
            #if 0.6 <= np.sum(avk) <= 1.4:
            counter=counter+1
            if np.mod(counter,8)==0:
                axs[1].plot(avk, o3_z/1e3, color=cmap(color_count*0.25+0.01))#label='z = '+f'{o3_z.sel(o3_p=avk.o3_p).values/1e3:.0f}'+' km'
                color_count = color_count +1
            else:
                if counter==1:
                    axs[1].plot(avk,o3_z/1e3, color='silver', label='AVKs')
                else:
                    axs[1].plot(avk, o3_z/1e3, color='silver')
        counter=0
        color_count2 = 0
        for j, avk_ndacc in enumerate(new['O3.MIXING.RATIO.VOLUME_EMISSION_AVK'].data):
            #if 0.6 <= np.sum(avk) <= 1.4:
            counter=counter+1
            if np.mod(counter,8)==0:
                axs[2].plot(avk_ndacc, new_alt/1e3, color=cmap(color_count2*0.25+0.01))#label='z = '+f'{o3_z.sel(o3_p=avk.o3_p).values/1e3:.0f}'+' km'
                color_count2 = color_count2 +1
            else:
                if counter==1:
                    axs[2].plot(avk_ndacc,new_alt/1e3, color='silver', label='AVKs')
                else:
                    axs[2].plot(avk_ndacc, new_alt/1e3, color='silver')
        
        axs[2].plot(0.0025*(100-new['O3.MIXING.RATIO.VOLUME_EMISSION_APRIORI.CONTRIBUTION'].data), new_alt/1e3, color='black')
        #axs[1].set_ylim(5,85)
        axs[1].set_xlim(-0.1,0.35)
        axs[1].set_title('Original GROMORA')
        axs[2].set_xlim(-0.1,0.35)
        axs[2].set_title('NDACC')
        axs[1].set_ylabel(y_lab)

        for ax in axs:
            ax.grid()
            ax.grid(which='both',  axis='x', linewidth=0.5)
            ax.grid(which='both',  axis='y', linewidth=0.5)

        plt.show()
    return new_ds

def GEOMS_2_NDACC(new_ds, outfolder):
    '''
    This function tries to convert the xarray dataset to HDF4 GEOMS format.
    It is derived from a code snippet sent by Ian Boyd and initially written
    by Bavo Langerock for FTIR data. THe full package can be found on:

    https://github.com/NCAR/sfit-processing-environment/tree/v3.0

    Adapted for GROMORA by E.S. 08.2022
    '''
    filename = outfolder+new_ds.FILE_NAME
    outputfile=os.path.splitext(filename)[0]+'.hdf'

    # Create the file:
    hdf4id=SD(outputfile,SDC.WRITE|SDC.TRUNC|SDC.CREATE)

    def dtype524(dtype5):
        if type(dtype5)!=str: dtype5=str(dtype5)
        if '|S' in dtype5: dtype4=SDC.CHAR8
        elif dtype5=='int8': dtype4=SDC.INT8
        elif dtype5=='int16': dtype4=SDC.INT16
        elif dtype5=='int32': dtype4=SDC.INT32
        elif dtype5=='float32': dtype4=SDC.FLOAT32
        elif dtype5=='float64': dtype4=SDC.FLOAT64
        return dtype4

    #WRITE ALL FILE attributes
    for att5,v in new_ds.attrs.items():
        #print('Writing file attribute %s=%s'%(att5,v))
        att4=hdf4id.attr(str(att5))
        if isinstance(v, str): dtype5='|S'
        else: dtype5=v.dtype
        if str(dtype5)=='|S0': continue #an empty attribute
        if att5=='FILE_NAME': v=os.path.basename(outputfile)
        try: att4.set(dtype524(dtype5),v)
        except: continue
    
    #WRITE Coordinates + ATTRIBUTES
    for var,v in new_ds.coords.items():
        #print('Writing variable %s with shape %s, dtype %s'%(var,v.shape,v.dtype))
        v=v.data[...] #get numpy arrays
        if v.shape==(): v=v.reshape(1,)
        vid4=hdf4id.create(str(var),dtype524(v.dtype),v.shape)
        if '|S' in str(v.dtype): v=np.array(v,dtype='S1')
        vid4[:]=v
        for att5,atv in new_ds[var].attrs.items():
            #print('\tSetting attribute %s (dtype=%s,shape=%s): %s'%(att5,atv.dtype,atv.shape,atv))
            att4=vid4.attr(str(att5))
            if isinstance(atv, str): dtype5='|S'
            else: dtype5=atv.dtype
            #dtype5=atv.dtype
            if 'float' in str(dtype5):
                try: atv=float(atv)
                except TypeError: atv=map(float,atv)
            if str(dtype5)=='|S0': continue #an empty attribute
            if 'int' in str(dtype5): atv = int(atv)
            att4.set(dtype524(dtype5),atv)
        vid4.endaccess()
    #WRITE VARIABLES + ATTRIBUTES
    for var,v in new_ds.items():
        #print('Writing variable %s with shape %s, dtype %s'%(var,v.shape,v.dtype))
        v=v.data[...] #get numpy arrays
        if v.shape==(): v=v.reshape(1,)
        vid4=hdf4id.create(str(var),dtype524(v.dtype),v.shape)
        if '|S' in str(v.dtype): v=np.array(v,dtype='S1')
        vid4[:]=v
        for att5,atv in new_ds[var].attrs.items():
            #print('\tSetting attribute %s (dtype=%s,shape=%s): %s'%(att5,atv.dtype,atv.shape,atv))
            att4=vid4.attr(str(att5))
            if isinstance(atv, str): dtype5='|S'
            else: dtype5=atv.dtype
            #dtype5=atv.dtype
            if 'float' in str(dtype5):
                try: atv=float(atv)
                except TypeError: atv=map(float,atv)
            if str(dtype5)=='|S0': continue #an empty attribute
            if 'int' in str(dtype5): atv = int(atv)
            att4.set(dtype524(dtype5),atv)
        vid4.endaccess()
    hdf4id.end()
    print('Written: '+outputfile)
    return

if __name__ == "__main__":

    write_new = True
    #instrument_name = 'GROMOS'
    d = pd.to_datetime(dt.datetime.now()-dt.timedelta(days=4))

    save_folder_gromos = '/storage/tub/instruments/gromos/NDACC/RapidDelivery/'
    save_folder_somora = '/storage/tub/instruments/somora/NDACC/'

    #dates = pd.date_range(start="2017-01-01",end="2017-01-01")

    if write_new:
        try:
            gromos = gromora_level2_GEOMS(instrument_name='GROMOS', date=d, spectros = ['AC240'] , ex = '_oper', new_z=1e3*np.arange(4, 92, 2), avk_corr=False, plot_tprofile=None, save_nc=False, RD=True)
            GEOMS_2_NDACC(gromos,outfolder=save_folder_gromos)
        except Exception as e:
            print(e)
        pass

    # plot_cycle = [12] # [8]
    # for d in dates:
    #     if write_new:
    #         try:
    #             gromos = gromora_level2_GEOMS(instrument_name='GROMOS', date=d, plot_tprofile=plot_cycle, save_nc=False)
    #             GEOMS_2_NDACC(gromos,outfolder=save_folder_gromos)
    #         except Exception as e:
    #             print(e)
    #             pass
    #         try:
    #             somora = gromora_level2_GEOMS(instrument_name='SOMORA', date=d, plot_tprofile=plot_cycle, save_nc=False)
    #             GEOMS_2_NDACC(somora,outfolder=save_folder_somora)
    #         except Exception as e:
    #             print(e)
    #             pass
    #     else:
    #         filename = 'groundbased_mwr.o3_ubern001_bern_20160101t000004z_20160101t235953z_012.nc'
    #         new_ds = xr.open_dataset(filename, engine='netcdf4')
    #         GEOMS_2_NDACC(new_ds,outfolder=save_folder_gromos)
