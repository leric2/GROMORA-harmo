"""
Batch procesing for extracting ECMWF data


"""
import os
import sys
sys.path.insert(0, '/home/es19m597/Documents/GROMORA/GROMORA-harmo/scripts/pyretrievals/')
#from retrievals.data.ecmwf import grib as ecmwf_grib
from retrievals.data.ecmwf import levels as ecmwf_levels
from glob import glob
import xarray as xr
import json
import logging
import numpy as np
import pandas as pd
from datetime import datetime

from ecmwf_extract_locations_v2 import extract_from_nc, find_files


if __name__ == '__main__':
    dateRange = pd.date_range(start='2022-12-01', end='2022-12-31')
    output_folder = '/storage/atmosphere/atmosphere/ecmwf/locations/'
    prefix = '/ecmwf_era5_'
    loc_file = '/storage/atmosphere/atmosphere/ecmwf/locations/locations_GROMORA_era5.json'

    fileList = []
    for d in dateRange:
        number_of_ECMWF_levels=137

        date = datetime.strftime(d, '%Y%m%d')
        basepath_location = os.path.join('/storage/atmosphere/atmosphere/ecmwf/era5/europe/',date[0:4])
        
        filename = basepath_location+'/ECMWF_ERA5_'+date+'.nc'
        fileList.append(filename)
    
    with open(loc_file) as f:
        locations = json.load(f)

    extract_from_nc(fileList, locations, output_folder, prefix, number_of_ECMWF_levels)