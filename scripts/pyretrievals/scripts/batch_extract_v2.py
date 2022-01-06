"""
Batch procesing for extracting ECMWF data


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
from datetime import datetime

from ecmwf_extract_locations_v2 import extract_from_nc, find_files


if __name__ == '__main__':
    dateRange = pd.date_range(start='2021-01-01', end='2021-11-25')
    output = '/storage/tub/instruments/gromos/ECMWF_Bern/ecmwf_oper_v2_'
    loc_file = '/home/esauvageat/Documents/ECMWF/locations_GROSOM.json'

    if dateRange[0]> datetime(2013,6,24):
        number_of_ECMWF_levels=137
    else:
        number_of_ECMWF_levels=91

    fileList = []
    for d in dateRange:
        date = datetime.strftime(d, '%Y%m%d')
        basepath_location = os.path.join('/storage/tub/atmosphere/ecmwf/oper/',date[0:4])
        
        filename = basepath_location+'/ECMWF_OPER_v2_'+date+'.nc'
        fileList.append(filename)
    
    with open(loc_file) as f:
        locations = json.load(f)

    extract_from_nc(fileList, locations, output, number_of_ECMWF_levels)

