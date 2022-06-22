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
import pandas as pd
import xarray as xr
from dotenv import load_dotenv
from pytz import timezone
from gromora_utils import *
from gromora_time import *

# # For ARTS, we need to specify some paths
# load_dotenv('/opt/anaconda/.env.birg-arts24')
# load_dotenv('/home/esauvageat/Documents/ARTS/.env.moench-arts2.4')

ARTS_DATA_PATH = os.environ['ARTS_DATA_PATH']
ARTS_BUILD_PATH = os.environ['ARTS_BUILD_PATH']
ARTS_INCLUDE_PATH = os.environ['ARTS_INCLUDE_PATH']

from retrievals import arts
from retrievals.data.ecmwf import ECMWFLocationFileStore
from retrievals.data.ecmwf import levels
from retrievals.data import interpolate
from retrievals.data import p_interpolate
from retrievals.arts.atmosphere import p2z_simple, z2p_simple
#from typhon.arts.xml import load
from pyarts.xml import load
from GROMORA_library import cmap

def gromora_level2_ndacc(instrument_name= "GROMOS", date= dt.date(2021, 6 , 27), plot_tprofile=False):
    """Function to convert GROMORA level 2 to GEOMS compliant NDACC data

    For now, it saves the new file as HDF5 file.

    TODO: change outputs to HDF4 !

    Args:
        instrument_name (str, optional): Name of the instrument. Defaults to "GROMOS".
        date (datetime, optional): the date to treat. Defaults to dt.date(2021, 6 , 27).
        plot_tprofile (bool, optional): a boolean to plot the comparisons between the og and the NDACC file. Defaults to False.

    Returns:
        xarray dataset: the new generated dataset, GEMOS compliant
    """
    #instrument_name = "GROMOS"

    #date = pd.date_range(start=sys.argv[1], end=sys.argv[2])
    #date = datetime.date(2016,1,2)
    #date = pd.date_range(start='2018-08-01', end='2018-08-02')
    #date = dt.date(2021, 6 , 27)
    #[pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=7)), pd.to_datetime(datetime.datetime.now()-datetime.timedelta(days=6))]

    integration_strategy = 'classic'
    spectros = ['AC240'] 
    int_time = 1
    ex = '_v2'
    # ex = '_waccm_low_alt'

    new_L2 = True

    #plotfolder = '/scratch/GROSOM/Level2/GROMORA_retrievals_v2/'
    plotfolder = '/storage/tub/instruments/gromos/level2/GROMORA/oper/'
    outfolder = '/home/es19m597/Documents/GROMORA/NDACC/'
    cont_name = 'h2o_continuum_x' 

    colormap = 'cividis'  # 'viridis' #, batlow_map cmap_crameri cividis

    measurement_response_limit = 0.75
    fill_value = -90000

    #############################################################################
    #############################################################################
    if instrument_name == "GROMOS":
        import gromos_classes as gromos
        basename_lvl1 = os.path.join(
            '/storage/tub/instruments/gromos/level1/GROMORA/v2/', str(date.year))
        #basename_lvl2 = "/scratch/GROSOM/Level2/GROMORA_retrievals_polyfit2/"
        if new_L2:
            basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v2/"+str(date.year)
        else:
            basename_lvl2 = "/storage/tub/instruments/gromos/level2/GROMORA/v1/"+str(date.year)
        instrument = gromos.GROMOS_LvL2(
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

    level2_dataset = instrument.read_level2(
        spectrometers=spectros,
        extra_base=ex
    )

    dataset = level2_dataset['AC240'].drop_dims(['f', 'h2o_continuum_p', 'h2o_continuum_p_avk', 'poly_order', 'f_shift_grid', 'sine_grid'])
    #dataset['o3_p'] = dataset.o3_z.mean(dim='time').data
    new_z = dataset.o3_z.mean(dim='time').data

    #interp_ds = dataset.interp(o3_p=z2p_simple(new_z))
    julian_dates = mjd2k_date(pd.to_datetime(dataset.time.data))

    #.replace(tzinfo=gromora_tz)

    lat = dataset.lat.mean(dim='time').data
    lon = dataset.lon.mean(dim='time').data  
    
    if not 'tropospheric_opacity' in list(dataset.data_vars.keys()):
        integrated_dataset, flags, integrated_meteo = instrument.read_level1b(
            no_flag=False, meta_data=True, extra_base='')
        level1b = integrated_dataset['AC240']
        altitude_instrument = np.median(level1b.alt.data)
    else:
        NotImplementedError('read tropospheric opactiy from level 1b !')

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
        lst, ha, sza, night, tc = get_LST_from_GROMORA(datetime64_2_datetime(t).replace(tzinfo=gromora_tz),  lat, lon, check_format=False)
        mean_sza.append(sza)
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
        #ds['temperature'] = (['time', 'o3_p'], ptz.t.data)

        # Interpolation of the dataset on the new altitude coordinates:
        t_interp = p_interpolate(ds.o3_p.data, ptz.p.values, ptz.t.values)
        ds.coords['alt'] = ('o3_p',ds.o3_z.data)
        ds = ds.swap_dims({'o3_p':'alt'})

        interp_ds = ds.interp(alt=new_z)

        # Clean ozone
        ozone = interp_ds.o3_x.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data
        N_o3 = o3_vmr2number_density(ozone, interp_ds.o3_p,t_interp)
        o3_number_density.append(N_o3)
        tco=0
        for i in range(len(N_o3)-1):
            tco = tco + 0.5*(np.nansum([N_o3[i],N_o3[i+1]]))*(new_z[i+1]-new_z[i])
        o3_total_column.append(tco/DU)

        o3.append(1e6*ozone)
        o3_p.append(1e2*interp_ds.o3_p)
        o3_random.append(1e6*interp_ds.o3_eo.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data)
        o3_sytematic.append(1e6*interp_ds.o3_es.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data)
        o3_tot.append(1e6*(interp_ds.o3_es+ interp_ds.o3_eo).where(interp_ds.o3_mr>measurement_response_limit, np.nan).data)
        o3_resolution.append(interp_ds.o3_fwhm.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data)
        o3_apriori.append(1e6*interp_ds.o3_xa)
        o3_apriori_contribution.append(100*(1-interp_ds.o3_mr.where(interp_ds.o3_mr>measurement_response_limit, np.nan).data))
        o3_avk.append(interp_ds.o3_avkm.data)
        temperature.append(t_interp)
        o3_mmr = o3_vmr2mmr(ozone)

        #tco = np.nansum(0.5*np.ediff1d(N_o3)*np.ediff1d(new_z))
        
    #o3_density = interp_ds.o3_p.data/t_interp/K_B*o3
    #o3_density = interp_ds.o3_p.data/t_interp/K_B*o3
    first_time = level1b.first_sky_time.reindex_like(dataset.time, method='nearest', tolerance='1H')
    last_time = level1b.last_sky_time.reindex_like(dataset.time, method='nearest', tolerance='1H')
    start_date_iso = pd.to_datetime(
        first_time[0].data
        ).round('1s').strftime('%Y%m%dT%H%M%Sz')

    stop_date_iso = pd.to_datetime(
        last_time[-1].data
        ).round('1s').strftime('%Y%m%dT%H%M%Sz')
    file_version='012'
    filename = 'groundbased_mwr.o3_'+instrument.affiliation+'_'+instrument.location.lower()+'_'+start_date_iso+'_'+stop_date_iso+'_'+file_version+'.hdf'
    
    # Global attributes dictionary:
    global_attrs = instrument.global_attributes_ndacc
    global_attrs['DATA_DISCIPLINE'] = 'ATMOSPHERIC.PHYSICS;REMOTE.SENSING;GROUNDBASED',
    global_attrs['DATA_GROUP'] = 'EXPERIMENTAL;PROFILE.STATIONARY',
    global_attrs['DATA_LOCATION'] = instrument.location,
    global_attrs['DATA_SOURCE']=instrument.source,
    global_attrs['DATA_START_DATE'] = start_date_iso
    global_attrs['DATA_STOP_DATE'] = stop_date_iso

    global_attrs['DATA_FILE_VERSION'] = file_version
    global_attrs['FILE_NAME']=filename,
    global_attrs['FILE_GENERATION_DATE'] = dt.datetime.now().isoformat(timespec='seconds')
    global_attrs['DATA_MODIFICATIONS'] = 'New harmonized retrievals from Swiss MWRs for FFTS time period',
    global_attrs['DATA_CAVEATS'] = 'Ozone profiles are estimated with the Optimal Estimation Method as implemented in the Atmospheric Radiative Transfer Simulator (ARTS)',
    global_attrs['DATA_QUALITY'] = '',
    global_attrs['DATA_TEMPLATE'] = 'GEOMS-TE-MWR-003',
    global_attrs['DATA_PROCESSOR'] = 'Harmonized time series processed with '+dataset.arts_version,
    global_attrs['FILE_ACCESS']='AVDC;NDACC',
    global_attrs['FILE_DOI']='',
    global_attrs['FILE_ASSOCIATION']='NDACC',
    global_attrs['FILE_META_VERSION']=''

    # var_dict = dict()
    # var_dict['LATITUDE.INSTRUMENT'] = (interp_ds.lat.data)
    # var_dict['LONGITUDE.INSTRUMENT'] = (interp_ds.lon.data),
    # var_dict['ALTITUDE.INSTRUMENT'] = (altitude_instrument),
    # var_dict['ANGLE.VIEW_AZIMUTH'] = ('DATETIME', dataset.obs_aa.data),
    # var_dict['ANGLE.VIEW_ZENITH'] = ('DATETIME', dataset.obs_za.data),
    # var_dict['ANGLE.SOLAR_ZENITH_MEAN'] = ('DATETIME',mean_sza),
    # var_dict['OPACITY.ATMOSPHERIC_EMISSION'] = ('DATETIME',level1b.tropospheric_opacity.reindex_like(dataset.time, method='nearest', tolerance='1H').data),
    # var_dict['DATETIME.START'] = ('DATETIME', mjd2k_date(pd.to_datetime(first_time.data))),
    # var_dict['DATETIME.STOP'] = ('DATETIME', mjd2k_date(pd.to_datetime(last_time.data))),
    # var_dict['INTEGRATION.TIME'] = ('DATETIME',instrument.cycle_duration*level1b.number_of_sky_spectra.reindex_like(dataset.time, method='nearest', tolerance='1H').data),
    # var_dict['PRESSURE_INDEPENDENT'] = (['DATETIME', 'ALTITUDE'], o3_p),
    # var_dict['TEMPERATURE_INDEPENDENT'] = (['DATETIME', 'ALTITUDE'], temperature),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.RANDOM.STANDARD'] = (['DATETIME', 'ALTITUDE'], o3_random ),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.SYTEMATIC.STANDARD']= (['DATETIME', 'ALTITUDE'], o3_sytematic ),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.COMBINED.STANDARD'] = (['DATETIME', 'ALTITUDE'], o3_tot ),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION'] = (['DATETIME', 'ALTITUDE'], o3 ),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION_RESOLUTION.ALTITUDE'] = (['DATETIME', 'ALTITUDE'], o3_resolution ),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION_APRIORI'] = (['DATETIME', 'ALTITUDE'], o3_apriori ),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION_APRIORI.CONTRIBUTION'] = (['DATETIME', 'ALTITUDE'], o3_apriori_contribution ),
    # var_dict['O3.MIXING.RATIO.VOLUME_EMISSION_AVK'] = (['DATETIME', 'ALTITUDE', 'ALTITUDE'], o3_avk ),
    # var_dict['O3.COLUMN.PARTIAL_EMISSION'] = (['DATETIME'], o3_total_column ),
    # var_dict['O3.NUMBER.DENSITY_EMISSION'] = (['DATETIME', 'ALTITUDE'], o3_apriori ),
    var_dict= dict(
            LATITUDE = (interp_ds.lat.data),
            LONGITUDE = (interp_ds.lon.data),
            ALTITUDE_INSTRUMENT = (altitude_instrument),
            ANGLE_VIEW_AZIMUTH = ('DATETIME', dataset.obs_aa.data),
            ANGLE_VIEW_ZENITH = ('DATETIME', dataset.obs_za.data),
            ANGLE_SOLAR_ZENITH_MEAN = ('DATETIME',mean_sza),
            #NOISE= ('DATETIME', dataset.median_noise.data),
            OPACITY_ATMOSPHERIC_EMISSION = ('DATETIME',level1b.tropospheric_opacity.reindex_like(dataset.time, method='nearest', tolerance='1H').data),
            DATETIME_START = ('DATETIME', mjd2k_date(pd.to_datetime(first_time.data))),
            DATETIME_STOP = ('DATETIME', mjd2k_date(pd.to_datetime(last_time.data))),
            INTEGRATION_TIME = ('DATETIME',instrument.cycle_duration*level1b.number_of_sky_spectra.reindex_like(dataset.time, method='nearest', tolerance='1H').data),
            PRESSURE_INDEPENDENT = (['DATETIME', 'ALTITUDE'], o3_p),
            TEMPERATURE_INDEPENDENT = (['DATETIME', 'ALTITUDE'], temperature),
            O3_MIXING_RATIO_VOLUME_EMISSION = (['DATETIME', 'ALTITUDE'], o3 ),
            O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_RANDOM_STANDARD = (['DATETIME', 'ALTITUDE'], o3_random ),
            O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_SYTEMATIC_STANDARD= (['DATETIME', 'ALTITUDE'], o3_sytematic ),
            O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_COMBINED_STANDARD = (['DATETIME', 'ALTITUDE'], o3_tot ),
            O3_MIXING_RATIO_VOLUME_EMISSION_RESOLUTION_ALTITUDE = (['DATETIME', 'ALTITUDE'], o3_resolution ),
            O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI = (['DATETIME', 'ALTITUDE'], o3_apriori ),
            O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI_CONTRIBUTION = (['DATETIME', 'ALTITUDE'], o3_apriori_contribution ),
            O3_MIXING_RATIO_VOLUME_EMISSION_AVK = (['DATETIME', 'ALTITUDE', 'ALTITUDE'], o3_avk ),
            O3_COLUMN_PARTIAL_EMISSION = (['DATETIME'], o3_total_column ),
            O3_NUMBER_DENSITY_EMISSION = (['DATETIME', 'ALTITUDE'], o3_apriori ),
        )
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

    # for s in list(new_ds.coords):
    #     n = list_indices[s]
    #     new_ds[s].attrs['VAR_NAME'] = s
    #     new_ds[s].attrs['VAR_DESCRIPTION'] = VAR_DESCRIPTION_List[n]
    #     new_ds[s].attrs['VAR_NOTES'] = VAR_NOTES_List[n]
    #     #new_ds[s].attrs['VAR_SIZE'] = VAR
    #     new_ds[s].attrs['VAR_DEPEND'] = VAR_DEPEND_List[n]
    #     new_ds[s].attrs['VAR_UNITS'] = VAR_UNITS_List[n]
    #     new_ds[s].attrs['VAR_SI_CONVERSION'] = VAR_SI_CONVERSION_List[n]
    #     new_ds[s].attrs['VAR_VALID_MIN'] = VAR_VALID_MIN_List[n]
    #     new_ds[s].attrs['VAR_VALID_MAX'] = VAR_VALID_MAX_List[n]
    #     new_ds[s].attrs['VAR_FILL_VALUE'] = -900000

    for s in list(new_ds.coords)+varname:
        n = list_indices[s]
        valid_min = VAR_VALID_MIN_List[n]
        valid_max = VAR_VALID_MAX_List[n]

        if s =='O3_NUMBER_DENSITY_EMISSION':
            fillV = -9e19
        else:
            fillV = -900000

        # Fill the nan values:
        if s != 'O3_MIXING_RATIO_VOLUME_EMISSION_AVK':
            new_ds[s].where(new_ds[s].data < valid_min, np.nan)
            new_ds[s].where(new_ds[s].data > valid_max, np.nan)
            new_ds[s] = new_ds[s].fillna(fillV)
        else:
            avk = new_ds[s].data
            avk[np.isnan(new_ds[s].data)] = fillV
            avk[avk<valid_min] = fillV
            avk[avk>valid_max] = fillV

            new_ds[s] = (['DATETIME', 'ALTITUDE', 'ALTITUDE'], avk)
        
        new_ds[s].attrs['VAR_NAME'] = SDS_Name_List[n]
        new_ds[s].attrs['VAR_DESCRIPTION'] = VAR_DESCRIPTION_List[n]
        new_ds[s].attrs['VAR_NOTES'] = VAR_NOTES_List[n]
        #new_ds[s].attrs['VAR_SIZE'] = VAR
        new_ds[s].attrs['VAR_DEPEND'] = VAR_DEPEND_List[n]
        new_ds[s].attrs['VAR_UNITS'] = VAR_UNITS_List[n]
        new_ds[s].attrs['VAR_SI_CONVERSION'] = VAR_SI_CONVERSION_List[n]
        new_ds[s].attrs['VAR_DATA_TYPE'] = 'DOUBLE'
        new_ds[s].attrs['VAR_VALID_MIN'] = valid_min
        new_ds[s].attrs['VAR_VALID_MAX'] = valid_max
        new_ds[s].attrs['VAR_FILL_VALUE'] = fillV
        new_ds[s].attrs['_FillValue'] = fillV
            
        if VAR_DEPEND_List[n] == 'CONSTANT':
            new_ds[s].attrs['VAR_SIZE'] = '1'
        elif VAR_DEPEND_List[n] == 'DATETIME':
            new_ds[s].attrs['VAR_SIZE'] = str(len(new_ds.DATETIME))
        elif VAR_DEPEND_List[n] == 'DATETIME;ALTITUDE':
            new_ds[s].attrs['VAR_SIZE'] = str(len(new_ds.DATETIME))+';'+str(len(new_ds.ALTITUDE))
        elif VAR_DEPEND_List[n] == 'DATETIME;ALTITUDE;ALTITUDE':
            new_ds[s].attrs['VAR_SIZE'] = str(len(new_ds.DATETIME))+';'+str(len(new_ds.ALTITUDE))+';'+str(len(new_ds.ALTITUDE))

        # rename var:
        new_ds = new_ds.rename_vars({s:SDS_Name_List[n]})
       
        variables_all = variables_all+SDS_Name_List[n]+';'

    new_ds.attrs['DATA_VARIABLE'] = variables_all

    new_ds.to_netcdf(outfolder+filename)

    if plot_tprofile is not None:
        fig, axs = plt.subplots(nrows=1, ncols=4, sharey=False, figsize=(24,16))
        i = plot_tprofile[0]

        ds = level2_dataset['AC240'].isel(time=i)
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

        #p_grid = np.logspace(-3, 5, 47)
        z_grid = p2z_simple(ds.o3_p.data)
        p_grid = z2p_simple(np.arange(1e3, 95e3, 2e3))
        p_grid = ds.o3_p.data

        z_interp = p_interpolate(p_grid, ptz.p.values, ptz.z.values )
        t_interp = p_interpolate(p_grid, ptz.p.values, ptz.t.values)
        axs[0].plot(ptz.t.values, ptz.z.values/1e3, 'x')
        axs[0].plot(t_interp, z_interp/1e3, 'o')
        axs[0].plot(new_ds.isel(DATETIME=i).TEMPERATURE_INDEPENDENT, new_ds.isel(DATETIME=i).ALTITUDE/1e3)
        axs[0].set_xlim(100,340)

        o3_og = ds.o3_x #.where(ds.o3_mr>measurement_response_limit, np.nan)
        o3 = new['O3.MIXING.RATIO.VOLUME_EMISSION']
        o3_z = ds.o3_z 
        new_alt = new.ALTITUDE
        o3_p = ds.o3_p
        y_lab = 'Altitude [km]'

        axs[3].plot(1e6*o3_og, o3_z/1e3, label='OG')
        axs[3].plot(o3, new_alt/1e3, label='NDACC')
        axs[3].legend()
        axs[3].set_xlim(-0.1,10)
        #axs[1].fill_betweenx(y_axis, (o3-error)*to_ppm,(o3+error)*to_ppm, color=col, alpha=0.5)
        #axs[1].plot(o3*to_ppm, y_axis,'-', linewidth=1.5, label='retrieved',color=col)
        #axs[0].plot(o3_apriori*to_ppm, y_axis, '--', linewidth=1.5, label='apriori',color='k')
        #axs[0].set_title('O$_3$ VMR')
        #axs[0].set_xlim(-0.5,9)

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

        #axs[1].set_ylim(5,85)
        axs[1].set_xlim(-0.1,0.3)
        axs[1].set_title('Original GROMORA')
        axs[2].set_xlim(-0.1,0.3)
        axs[2].set_title('NDACC')
        axs[1].set_ylabel(y_lab)

        for ax in axs:
            ax.grid()
            ax.grid(which='both',  axis='x', linewidth=0.5)
            ax.grid(which='both',  axis='y', linewidth=0.5)

        plt.show()
    return new_ds

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

SDS_Name_List = ['LATITUDE.INSTRUMENT',
                 'LONGITUDE.INSTRUMENT',
                 'ALTITUDE.INSTRUMENT',
                 'DATETIME',
                 'ANGLE.VIEW_AZIMUTH',
                 'ANGLE.VIEW_ZENITH_MEAN',
                 'ANGLE.SOLAR_ZENITH_MEAN',
                 'OPACITY.ATMOSPHERIC_EMISSION',
                 'DATETIME.START',
                 'DATETIME.STOP',
                 'INTEGRATION.TIME',
                 'ALTITUDE',
                 'PRESSURE_INDEPENDENT',
                 'TEMPERATURE_INDEPENDENT',
                 'O3.MIXING.RATIO.VOLUME_EMISSION',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.RANDOM.STANDARD',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.SYSTEMATIC.STANDARD',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.COMBINED.STANDARD',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_RESOLUTION.ALTITUDE',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_APRIORI',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_APRIORI.CONTRIBUTION',
                 'O3.MIXING.RATIO.VOLUME_EMISSION_AVK',
                 'O3.COLUMN.PARTIAL_EMISSION',
                 'O3.NUMBER.DENSITY_EMISSION']

SDS_DataType_List = ['float',
                     'float',
                     'int32',
                     'double',
                     'float',
                     'float',
                     'float',
                     'float',
                     'double',
                     'double',
                     'double',
                     'int32',
                     'float',
                     'float',
                     'float',
                     'float',
                     'float',
                     'float',
                     'int32',
                     'float',
                     'float',
                     'float',
                     'float',
                     'float']

VAR_DESCRIPTION_List = ['Latitude of the observation site (deg)',
                        'Longitude of the observation site (deg)',
                        'Altitude of the observation site (m)',
                        'Representative Date and Time of the Measurement (MJD2K)',
                        'Viewing Azimuth angle (deg)',
                        'Mean Viewing Zenith angle (deg)',
                        'Mean Solar Zenith angle of measurement (deg)',
                        'Tropospheric Optical Depth at 142 GHz (Np)',
                        'Start Date and Time of the Measurement (MJD2K)',
                        'Stop Date and Time of the Measurement (MJD2K)',
                        'Effective Measurement Time (h)',
                        'Altitude as geometric height (m)',
                        'Derived Atmospheric Pressure (hPa)',
                        'Derived Temperature Profile (K)',
                        'O3 Volume Mixing Ratio Profile (ppmv)',
                        'O3 Profile Observation Error (ppmv)',
                        'O3 Profile Smoothing Error (ppmv)',
                        'O3 Profile Total Error (ppmv)',
                        'O3 Mixing Ratio Altitude Resolution (m)',
                        'A priori O3 Profile (ppmv)',
                        'Contribution from the A priori to the retrieved O3 Profile (%)',
                        'O3 Retrieval Averaging Kernels',
                        'O3 partial column density (DU)',
                        'O3 number density (molecules m-3)'];
                        
VAR_NOTES_List = [' ',
                  ' ',
                  ' ',
                  'Averaged time of sky measurements',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  'Merge of ECMWF operational analysis and CIRA86',
                  'Merge of ECMWF and CIRA86',
                  ' ',
                  ' ',
                  ' ',
                  ' ',
                  'Resolution defined as FWHM of the averaging kernel peak for that altitude level.',
                  'Monthly day and night profile extracted from WACCM free running model',
                  ' ',
                  ' ',
                  ' ',
                  ' ']
                    
VAR_DEPEND_List = ['CONSTANT',             
                   'CONSTANT',
                   'ALTITUDE',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'DATETIME',
                   'ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE',
                   'DATETIME;ALTITUDE;ALTITUDE',
                   'DATETIME',
                   'DATETIME;ALTITUDE']
               
VAR_UNITS_List = ['deg',
                  'deg',
                  'm',
                  'MJD2K',
                  'deg',
                  'deg',
                  'deg',
                  'Np',
                  'MJD2K',
                  'MJD2K',
                  'h',
                  'm',
                  'hPa',
                  'K',
                  'ppmv',
                  'ppmv',
                  'ppmv',
                  'ppmv',
                  'm',
                  'ppmv',
                  '%',
                  '1',
                  'DU',
                  'molec m-3']
    
VAR_SI_CONVERSION_List = ['0.0;1.74533E-2;rad',
                          '0.0;1.74533E-2;rad',
                          '0;1;m',
                          '0.0;86400.0;s',
                          '0.0;1.74533E-2;rad',
                          '0.0;1.74533E-2;rad',
                          '0.0;1.74533E-2;rad',
                          '0.0;1.0;1',
                          '0.0;86400.0;s',
                          '0.0;86400.0;s',
                          '0.0;3600.0;s',
                          '0;1;m',
                          '0.0;1.0E2;kg m-1 s-2',
                          '0.0;1.0;K',
                          '0.0;1.0E-6;1',
                          '0.0;1.0E-6;1',
                          '0.0;1.0E-6;1',
                          '0.0;1.0E-6;1',
                          '0;1;m',
                          '0.0;1.0E-6;1',
                          '0.0;0.01;1',
                          '0.0;1.0;1',
                          '0.0;4.4614E-4;mol m-2',
                          '0.0;1.66054E-24;mol m-3']

VAR_VALID_MIN_List = [-90,
                      -180,
                      0,
                      -2200,
                      0,
                      0,
                      0,
                      0,
                      np.nan,
                      np.nan,
                      0,
                      0,
                      0,
                      100,
                      -0.1,
                      0.1,
                      0.1,
                      0.1,
                      5000,
                      -0.1,
                      -70,
                      -1,
                      0,
                      -9.e19]

VAR_VALID_MAX_List = [90,
                      180,
                      5000,
                      np.nan,
                      360,
                      180,
                      180,
                      5,
                      np.nan,
                      np.nan,
                      1000,
                      120000,
                      500,
                      500,
                      20,
                      20,
                      20,
                      20,
                      35000,
                      20,
                      200,
                      1,
                      500,
                      1.e19]

list_indices = {
    'DATETIME':3, 
    'ALTITUDE':11, 
    'LATITUDE':0,
    'LONGITUDE':1, 
    'ALTITUDE_INSTRUMENT':2,
    'ANGLE_VIEW_AZIMUTH':4,
    'ANGLE_VIEW_ZENITH':5,
    'ANGLE_SOLAR_ZENITH_MEAN':6,
    'OPACITY_ATMOSPHERIC_EMISSION':7,
    'DATETIME_START':8,
    'DATETIME_STOP':9,
    'INTEGRATION_TIME':10,
    'PRESSURE_INDEPENDENT':12,
    'TEMPERATURE_INDEPENDENT':13,
    'O3_MIXING_RATIO_VOLUME_EMISSION':14,
    'O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_RANDOM_STANDARD':15,
    'O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_SYTEMATIC_STANDARD':16,
    'O3_MIXING_RATIO_VOLUME_EMISSION_UNCERTAINTY_COMBINED_STANDARD':17,
    'O3_MIXING_RATIO_VOLUME_EMISSION_RESOLUTION_ALTITUDE':18,
    'O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI':19,
    'O3_MIXING_RATIO_VOLUME_EMISSION_APRIORI_CONTRIBUTION':20,
    'O3_MIXING_RATIO_VOLUME_EMISSION_AVK':21,
    'O3_COLUMN_PARTIAL_EMISSION':22,
    'O3_NUMBER_DENSITY_EMISSION':23
    }

if __name__ == "__main__":
    instrument_name = 'GROMOS'
    d = dt.date(2020, 1 , 27)
    
    dates = pd.date_range(start="2019-09-09",end="2019-09-09")

    plot_cycle = None # [8]
    for d in dates:
        new_ds = gromora_level2_ndacc(instrument_name=instrument_name, date=d, plot_tprofile=plot_cycle)
