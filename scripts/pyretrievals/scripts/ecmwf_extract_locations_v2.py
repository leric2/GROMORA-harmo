"""
Extract ecmwf OPER v2 data for specified locations from global netCDF or GRIB2 files.
An update from the OG function from Jonas for GRIB1 files

For list_params to work, the file has to be in the same directory...

Usage:
  ecmwf_extract_locations.py extract <locations_file> <basepath> <output_prefix> <file_type> <file>...
  ecmwf_extract_locations.py list_params <file_type> <file>

Options:
  -h --help     Show this screen.
  --version     Show version.
"""
import os
#from retrievals.data.ecmwf import grib as ecmwf_grib
from retrievals.data.ecmwf import levels as ecmwf_levels
from glob import glob
import xarray as xr
from docopt import docopt
import json
import logging
#import pygrib
import numpy as np
import pandas as pd


logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.DEBUG)


ATTRS = {
    'lat': {
        'standard_name': 'latitude',
        'long_name': 'Latitude',
        'units': 'degree_north',
        'axis': 'Y',
    },
    'lon': {
        'standard_name': 'longitude',
        'long_name': 'Longitude',
        'units': 'degree_east',
        'axis': 'X',
    },
    'time': {
        'standard_name': 'time',
        'long_name': 'Time',
    },
    'level': {
        'long_name': 'ECMWF model level'
    },
    'pressure': {
        'long_name': 'Pressure',
        'standard_name': 'air_pressure',
        'units': 'Pa',
        'axis': 'Z',
        'positive': 'down',
    }
}

def find_files(basepath, file):
    if len(file) > 1:
        return sorted(os.path.join(basepath, file))
    return sorted(glob(os.path.join(basepath, file[0])))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def extract_from_grib(grib_files, locations, output_prefix):
    parameters, dims = get_dims(grib_files[0])
    n_times = len(dims['time'])
    n_levels = len(dims['level'])

    coerced_coords = dict()
    data = dict()
    for loc_name, loc in locations.items():
        # find coerced coords
        lat = find_nearest(dims['lat'], loc['lat'])
        lon = find_nearest(dims['lon'], loc['lon'])
        i_lat = np.where(dims['lat'] == lat)
        i_lon = np.where(dims['lon'] == lon)
        coerced_coords[loc_name] = (lat, lon, i_lat, i_lon)

        # initialize arrays
        data[loc_name] = dict()
        for p in parameters:
            if p in ['25', 'Geopotential']:
                a = np.ndarray((1, n_times), dtype=np.float32)
            else:
                a = np.ndarray((n_levels, n_times), dtype=np.float32)
            data[loc_name][p] = a

    # Read every file
    for file in grib_files:
        logging.debug('Clear arrays')
        for d in data.values():
            for a in d.values():
                a.fill(np.nan)

        dims['time'] = []

        logging.info('Open ' + file)
        grbs = pygrib.open(file)
        for msg in grbs:
            if msg.analDate not in dims['time']:
                dims['time'].append(msg.analDate)
            i_time = dims['time'].index(msg.analDate)
            i_level = dims['level'].index(msg.level)
            values = msg.values

            # extract data for each location
            for loc_name, ps in data.items():
                lat, lon, i_lat, i_lon = coerced_coords[loc_name]
                v = values[i_lat, i_lon]
                ps[msg.parameterName][i_level, i_time] = v


        # Complie data and store
        for loc_name, xs in data.items():
            lat, lon, _, _ = coerced_coords[loc_name]
            data_vars = dict()
            for p, a in xs.items():
                slug = p.lower().replace(' ', '_')
                if p in ['Geopotential']:
                    data_vars[slug] = (('loc', 'time',),
                                       a[0, :][np.newaxis, :],
                                       {'grib_name': p})
                elif p in ['25']:
                    data_vars['logarithm_of_surface_pressure'] = (('loc', 'time',),
                                       a[0, :][np.newaxis, :],
                                       {'grib_name': p})
                else:
                    data_vars[slug] = (('loc', 'level', 'time'),
                                       a[np.newaxis, :, :],
                                       {'grib_name': p})
            # calculate pressure
            pressure = ecmwf_levels.pressure_levels(xs['25'])
            data_vars['pressure'] = (('loc', 'level', 'time'),
                                     pressure,
                                     ATTRS['pressure'])

            coords = {
                'time': dims['time'],
                'level': dims['level'],
                'lat': ('loc', [lat], ATTRS['lat']),
                'lon': ('loc', [lon], ATTRS['lon']),
                'loc': ('loc', [loc_name], {'long_name': 'Location identifier'}),
            }
            attrs = {
                'grib_file': file,
            }
            ds = xr.Dataset(data_vars, coords, attrs=attrs)
            ds['time'].encoding['units'] = 'seconds since 1970-01-01 00:00:00'

            day_str = pd.Timestamp(ds['time'].values[0]).strftime('%Y%m%d')
            out_fn = output_prefix + loc_name + '_' + day_str + '.nc'
            logging.info('Write ' + out_fn)
            ds.to_netcdf(out_fn, unlimited_dims=['time',])

def extract_from_nc(nc_files, locations, output_prefix, levels_n=137):
    '''
    Main function to extract ECMWF data at specified locations. Write 1 netCDF file per location.
    No interpolation done, only take the nearest existing datapoint using lat and lon.

    Args:
        nc_files (list): list of filename to read
        locations (dict): locations
        output_prefix (str): the prefix for the output netCDF file

    Returns:
        None
    '''
    for file in nc_files:
        ecmwf_nc = xr.open_dataset(
            file,
            mask_and_scale=True,
            decode_times=True,
            decode_coords=True,
            #use_cftime=True,
        )

        dims = ecmwf_nc.coords
    
        coerced_coords = dict()
        parameters = []
        for i in ecmwf_nc.data_vars.keys():
            #print(i)
            parameters.append(str(i))
            #p = ecmwf_nc[i].attrs['long_name']
            #slug = p.lower().replace(' ', '_')
            #parameters.append(p)

        for loc_name, loc in locations.items():
            # find coerced coords
            lat = find_nearest(dims['latitude'], loc['lat'])
            lon = find_nearest(dims['longitude'], loc['lon'])
            i_lat = np.where(dims['latitude'] == lat)
            i_lon = np.where(dims['longitude'] == lon)
            coerced_coords[loc_name] = (lat, lon, i_lat, i_lon)

            ds_loc = ecmwf_nc.sel(latitude = lat, longitude =lon)

            # Creating data_vars dictionnary:
            data_vars = dict()
            for p in parameters:
                var_name = ds_loc[p].attrs['long_name']
                slug = var_name.lower().replace('(','').replace(')','').replace(' ', '_')
                if p in ['lnsp','z']:
                    data_vars[slug] = (('loc','time'),ds_loc[p].values[:,0][np.newaxis, :],{'og_name': p})
                else:
                    #
                    data_vars[slug] = (
                        ('loc', 'level', 'time'),
                        np.transpose(ds_loc[p].values)[ np.newaxis, :, :],
                        {'og_name': p}
                    )

            # Calculate pressure
            if levels_n==137:
                pressure = ecmwf_levels.pressure_levels(ds_loc['lnsp'].values[:,0])
            elif levels_n==91:
                pressure = ecmwf_levels.pressure_levels_91(ds_loc['lnsp'].values[:,0])
            data_vars['pressure'] = (
                ('loc', 'level', 'time'),
                pressure[np.newaxis,:,:],
                ATTRS['pressure']
            )

            # Building a new xarray 
            coords = {
                'time': dims['time'],
                'level': dims['level'],
                'lat': ('loc', [lat], ATTRS['lat']),
                'lon': ('loc', [lon], ATTRS['lon']),
                'loc': ('loc', [loc_name], {'long_name': 'Location identifier'}),
            }
            attrs = {
                'nc_file': file,
            }

            ds = xr.Dataset(data_vars, coords, attrs=attrs)
            ds['time'].encoding['units'] = 'seconds since 1970-01-01 00:00:00'

            day_str = pd.Timestamp(ds['time'].values[0]).strftime('%Y%m%d')
            out_fn = output_prefix + loc_name + '_' + day_str + '.nc'
            logging.info('Write ' + out_fn)
            ds.to_netcdf(out_fn, unlimited_dims=['time',])

def get_dims(grib_file):
    ge = ecmwf_grib.GribECMWF(grib_file)

    msg = ge.grbs[1]
    lons = np.linspace(-float(msg['longitudeOfFirstGridPointInDegrees']),
                       float(msg['longitudeOfLastGridPointInDegrees']),
                       int(msg['Ni']))
    lats = np.linspace(float(msg['latitudeOfFirstGridPointInDegrees']),
                       float(msg['latitudeOfLastGridPointInDegrees']),
                       int(msg['Nj']))

    parameters = sorted(ge.index.sel_values('parameterName'))
    dims = {
        'level': sorted(ge.index.sel_values('level')),
        'time': sorted(ge.index.sel_values('analDate')),
        'lat': lats,
        'lon': lons,
    }

    return parameters, dims

def list_params(file_type, file):
    if file_type == 'GRIB2':
        logging.info('Create index for ' + str(file))
        ge = ecmwf_grib.GribECMWF(file)
        parameters = ge.index.sel_values('parameterName')
        print()
        for p in parameters:
            n_levels = len(ge.index.sel(parameterName=p)) / 4
            print(p, int(n_levels), 'levels')
    elif file_type == 'netCDF':
        ecmwf_nc = xr.open_dataset(
            file,
            mask_and_scale=True,
            decode_times=True,
            decode_coords=True,
            )
        for i in ecmwf_nc.data_vars.keys():
            print(ecmwf_nc[str(i)].attrs['long_name'], len(ecmwf_nc[i].level), 'levels')

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Extract ECMWF 2.0')

    if arguments['extract']:
        with open(arguments['<locations_file>']) as f:
            locations = json.load(f)
        if arguments['<file_type>'] == 'GRIB2':
            extract_from_grib(find_files(arguments['<basepath>'], arguments['<file>']),
                locations,
                arguments['<output_prefix>'])
        elif arguments['<file_type>'] == 'netCDF':
            extract_from_nc(find_files(arguments['<basepath>'], arguments['<file>']),
                locations,
                arguments['<output_prefix>'])
        else:
            raise ValueError('File type not implemented (GRIB2 or netCDF only)')    
    elif arguments['list_params']:
        list_params(arguments['<file_type>'], arguments['<file>'][0])
