#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 

@author: eric

Collection of functions for dealing with time


Including : 
    * a-priori data
"""
import os
import numpy as np
import retrievals
import xarray as xr
import pandas as pd
import math
import netCDF4
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from pytz import timezone, utc

from pysolar import solar

gromora_tz = timezone('Europe/Zurich')

def solar_zenith_angle(ha, lst, lat):
    # sun declination
    declination = -23.44*np.cos(np.deg2rad((pd.to_datetime(lst).dayofyear+10)*360/365))
    cos_declination = np.cos(np.deg2rad(declination))
    sin_declination = np.sin(np.deg2rad(declination))

    # Hour angle 
    cos_solar_hour_angle = np.cos(np.deg2rad(ha))

    # Sunrise/Sunset:
    cos_hour_angle_night = -np.tan(np.deg2rad(lat))*np.tan(np.deg2rad(declination))
    if cos_solar_hour_angle < cos_hour_angle_night:
        night = True
    else:
        night = False
    cos_sza = sin_declination*np.sin(np.deg2rad(lat)) + np.cos(np.deg2rad(lat))*cos_declination*cos_solar_hour_angle
    sza = np.rad2deg(np.arccos(cos_sza))
    return sza, night

def pysolar_sza(date, lat, lon):
    # Using pysolar package:
    time = utc.localize(datetime64_2_datetime(date))
    sza_pysolar = solar.get_altitude(lat, lon, time)
    print('solar elevation angle = ', sza_pysolar)
    return sza_pysolar

def datetime64_2_datetime(dt64):
    unix_epoch = np.datetime64(0, 's')
    one_second = np.timedelta64(1, 's')
    seconds_since_epoch = (dt64 - unix_epoch) / one_second
    dt = datetime.utcfromtimestamp(seconds_since_epoch)
    return dt

def equation_of_time(doy):
    B=(360/365)*(doy-81)
    eot = 9.87*np.sin(np.deg2rad(2*B)) - 7.53*np.cos(np.deg2rad(B)) -1.5*np.sin(np.deg2rad(B))
    return eot

def time_correction_factor(lon, lstm, eot):
    return 4*(lon - lstm) + eot

def local_solar_time(local_time, tc):
    return local_time + tc/60

def hour_angle(lst):
    hours_from_midnight = (lst - lst.replace(hour=0, minute=0,second=0, microsecond=0)).total_seconds()/3600
    return 15*(hours_from_midnight-12)

def get_LST_from_UTC(date, lat, lon):
    dt = utc.localize(datetime64_2_datetime(date))
    print('UTC time: ',dt) 
    local_time =  dt.astimezone(gromora_tz)
    print('Local time: ',local_time)

    doy = pd.to_datetime(dt).dayofyear

    eot = equation_of_time(doy)

    lstm = 15*local_time.utcoffset().seconds/3600
    tc = time_correction_factor(lon, lstm, eot)

    seconds_from_midnight = (local_time - local_time.replace(hour=0,minute=0,second=0, microsecond=0)).total_seconds()
    minutes_from_midnight = seconds_from_midnight / 60

    lst = local_time + timedelta(hours=tc/60)
    ha = hour_angle(lst)

    # print('Local solar time: ',lst)
    # print('Hour angle: ', str(ha))

    sza, night = solar_zenith_angle(ha, lst, lat)

    #print('solar elevation angle : ', 90-sza)

    return lst, ha, sza, night

if __name__ == "__main__":
    #date = spectro_dataset.time.data[0]
    date = np.datetime64('2017-12-30 08:12:20.123456')

    lat = np.array(46.95)
    lon = np.array(7.44)

    sza_pysolar = pysolar_sza(date, lat, lon)

    lst, ha, sza, night = get_LST_from_UTC(date, lat, lon)

    if night:
        print('Night !')
    else:
        print('Day !')

    
