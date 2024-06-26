"""
Batch procesing for extracting ECMWF data


"""
import sys, os
from os.path import dirname, abspath, join

sys.path.insert(0, '/home/es19m597/Documents/GROMORA/GROMORA-harmo/scripts/pyretrievals/')
sys.path.append(join(dirname(sys.path[0]),'pyretrievals'))

#from retrievals.data.ecmwf import grib as ecmwf_grib
from retrievals.data.ecmwf import levels as ecmwf_levels
from glob import glob
import xarray as xr
#from docopt import docopt
import json
import logging
#import pygrib
import numpy as np
import pandas as pd
from datetime import datetime

from ecmwf_extract_locations_v2 import extract_from_nc, find_files


if __name__ == '__main__':
    dateRange = pd.date_range(start='2023-04-10', end='2023-04-10')
    output = '/storage/atmosphere/atmosphere/ecmwf/locations/'
    loc_file = '/storage/atmosphere/atmosphere/ecmwf/locations/locations.json'
    prefix = '/ecmwf_oper_v2_'

    if dateRange[0]> datetime(2013,6,24):
        number_of_ECMWF_levels=137
    else:
        number_of_ECMWF_levels=91

    fileList = []
    for d in dateRange:
        date = datetime.strftime(d, '%Y%m%d')
        basepath_location = os.path.join('/storage/atmosphere/atmosphere/ecmwf/oper/',date[0:4])
        
        filename = basepath_location+'/ECMWF_OPER_v2_'+date+'.nc'
        fileList.append(filename)
    
    with open(loc_file) as f:
        locations = json.load(f)

    extract_from_nc(fileList, locations, output, prefix, number_of_ECMWF_levels)

